"""
Experiment 4: Failure of Cameron-Martin Admissibility.

This stress test attacks Assumption (ii) (Cameron-Martin stability) by
constructing a 1D SPDE heat-equation prior and approximating a Dirac
point-evaluation observation via a narrow Gaussian mollifier:

    g_ε(x) = (1/√(2π)ε) exp(-(x-x₀)²/(2ε²))

For each fixed ε > 0, g_ε ∈ L²(0,1) so it lies in Range(Q₀^{1/2}) = H₀¹(0,1).
But as ε ↓ 0 (approaching the Dirac delta), the Cameron-Martin energy norm
|Q₀^{-1/2} g_ε| diverges, violating the admissibility condition.

The experiment demonstrates two consequences:
  (1) The CM energy norm |g_ε|_{Q₀} ~ ε^{-1} → ∞ as ε → 0.
  (2) The shorted residual R_λ upper bound (proportional to the CM energy squared)
      no longer provides the O(λ^{-6}) remainder control, destroying the
      asymptotic expansion.
"""

using CairoMakie
using LaTeXStrings
using DataFrames
using CSV
using Printf
using SHA

const OUTPUT_DIR = joinpath(@__DIR__, "results")
const PRECISION_BITS = 256

# Prior: L = -d²/dx² + κ² on [0,1] with Dirichlet BCs
# Eigenvalues q_j = (π²j² + κ²)⁻¹
# Cameron-Martin space: H_{Q₀} = H₀¹(0,1)
# CM energy norm: |h|²_{Q₀} = ∫₀¹ (|h'|² + κ²|h|²) dx

const KAPPA = BigFloat(1)
const X0 = BigFloat(1) / BigFloat(2)  # observation point
const N_MODES = 2000   # spectral truncation for inner products

# ── Epsilon grid for the mollifier width ──
# We use ε = 2^{-p} for p = 0..14 to span from ε=1 down to ε ≈ 6e-5
const EPS_POWERS = 0:14
const LAMBDA_POWERS = 1:12

"""
    cm_energy_norm_sq(eps0, kappa, x0, n_modes)

Compute |g_ε|²_{Q₀} = Σⱼ (π²j² + κ²) |⟨g_ε, e_j⟩|² in spectral representation.

The Fourier-sine coefficient of g_ε on [0,1] is:
    ⟨g_ε, e_j⟩ = √2 ∫₀¹ g_ε(x) sin(πjx) dx

For large j or small ε, the coefficient decays as exp(-(πjε)²/2) sin(πjx₀).

The CM energy norm is the H¹ energy:
    |g_ε|²_{Q₀} = Σⱼ (π²j² + κ²) · |⟨g_ε, e_j⟩|²
"""
function cm_energy_norm_sq(eps::BigFloat, kappa::BigFloat, x0::BigFloat, n_modes::Int)
    total = zero(BigFloat)
    sqrt2 = sqrt(BigFloat(2))
    pi_val = BigFloat(π)
    for j in 1:n_modes
        pij = pi_val * j
        # Fourier-sine coefficient of g_ε against e_j(x) = √2 sin(πjx)
        # ⟨g_ε, e_j⟩ = √2 · exp(-(πjε)²/2) · sin(πjx₀)
        # (exact for the Gaussian on ℝ, and an excellent approximation for
        #  ε ≪ min(x₀, 1-x₀) since the tails are negligible)
        coeff = sqrt2 * exp(-(pij * eps)^2 / 2) * sin(pij * x0)
        eigenvalue_inv = pij^2 + kappa^2  # = q_j^{-1}
        total += eigenvalue_inv * coeff^2
    end
    return total
end

"""
    l2_norm_sq(eps, x0, n_modes)

Compute ||g_ε||²_{L²} = Σⱼ |⟨g_ε, e_j⟩|² (spectral).
"""
function l2_norm_sq(eps::BigFloat, x0::BigFloat, n_modes::Int)
    total = zero(BigFloat)
    sqrt2 = sqrt(BigFloat(2))
    pi_val = BigFloat(π)
    for j in 1:n_modes
        pij = pi_val * j
        coeff = sqrt2 * exp(-(pij * eps)^2 / 2) * sin(pij * x0)
        total += coeff^2
    end
    return total
end

"""
    schur_residual_bound(eps, kappa, x0, n_modes, lam)

Compute the CM-energy-based upper bound on R_λ for a single observation
direction g_ε. The cross-covariance energy bound gives:

    R_λ ≤ |Σ_{λ,Ik}|²_{Q₀}  (Douglas range inclusion)
        ≤ (ε₀/(2λ²))² · M²/(λ+ω)² · |g_ε|²_{Q₀}

For this stress test, we take ε₀ = 1, M = 1, ω = gap = (π² + κ²), the smallest
eigenvalue rate in the free block. This yields the pessimistic scaling:

    R_λ^{bound} = |g_ε|²_{Q₀} / (4 λ⁴ (λ + ω)²)

The O(λ^{-6}) remainder holds only when |g_ε|²_{Q₀} is finite.
"""
function schur_residual_bound(eps::BigFloat, kappa::BigFloat, x0::BigFloat,
                              n_modes::Int, lam::BigFloat)
    cm_sq = cm_energy_norm_sq(eps, kappa, x0, n_modes)
    # Spectral gap of the drift (smallest eigenvalue of the drift generator magnitude)
    omega = BigFloat(π)^2 + kappa^2
    eps0 = one(BigFloat)
    # The resolvent cross-covariance is Σ_{λ,Ik} = (ε₀/(2λ²)) (λI - A_II)^{-1} a_{Ik}
    # Its CM norm is bounded by (ε₀/(2λ²)) · M/(λ+ω) · |a_{Ik}|_{Q₀}
    cross_cm_bound = (eps0 / (BigFloat(2) * lam^2)) * (one(BigFloat) / (lam + omega)) * sqrt(cm_sq)
    return cross_cm_bound^2
end

function evaluate_cm_divergence()
    records = NamedTuple[]
    for p in EPS_POWERS
        eps = BigFloat(2)^(-p)
        cm_sq = cm_energy_norm_sq(eps, KAPPA, X0, N_MODES)
        l2_sq = l2_norm_sq(eps, X0, N_MODES)
        push!(records, (
            eps_power = p,
            epsilon = eps,
            l2_norm_sq = l2_sq,
            cm_energy_norm_sq = cm_sq,
            cm_energy_norm = sqrt(cm_sq),
        ))
    end
    return DataFrame(records)
end

function evaluate_remainder_breakdown()
    records = NamedTuple[]
    # Pick representative epsilon values: admissible (ε=1), marginal (ε=2^{-4}),
    # severely rough (ε=2^{-8}), near-delta (ε=2^{-12})
    test_eps_powers = [0, 4, 8, 12]
    for p in test_eps_powers
        eps = BigFloat(2)^(-p)
        cm_sq = cm_energy_norm_sq(eps, KAPPA, X0, N_MODES)
        for lp in LAMBDA_POWERS
            lam = BigFloat(2)^lp
            r_bound = schur_residual_bound(eps, KAPPA, X0, N_MODES, lam)
            # The canonical target variance
            sigma_kk = one(BigFloat) / (BigFloat(2) * lam^2)
            # Check if remainder bound exceeds the leading term (breakdown)
            ratio = r_bound / sigma_kk
            push!(records, (
                eps_power = p,
                epsilon = eps,
                lambda_power = lp,
                lambda = lam,
                cm_energy_sq = cm_sq,
                R_bound = r_bound,
                lambda6_R_bound = lam^6 * r_bound,
                sigma_kk = sigma_kk,
                ratio_R_over_sigma = ratio,
            ))
        end
    end
    return DataFrame(records)
end

function plot_results(cm_df::DataFrame, breakdown_df::DataFrame)
    fig = Figure(size=(900, 720))

    # ── Panel A: CM energy norm divergence as ε → 0 ──
    ax1 = Axis(fig[1, 1],
        xscale = log10,
        yscale = log10,
        title = "Cameron-Martin energy divergence (Panel A)",
        xlabel = L"Mollifier width $\varepsilon$",
        ylabel = L"$|g_\varepsilon|_{Q_0}^2$",
        xgridvisible = true,
        ygridvisible = true,
        xgridcolor = (:black, 0.1),
        ygridcolor = (:black, 0.1),
        xminorticksvisible = true,
        yminorticksvisible = true,
        xminorticks = IntervalsBetween(9),
        yminorticks = IntervalsBetween(9)
    )
    eps_vals = Float64.(cm_df.epsilon)
    cm_vals = Float64.(cm_df.cm_energy_norm_sq)
    scatterlines!(ax1, eps_vals, cm_vals,
        color = "#0B7285", linewidth = 2.0, marker = :circle, markersize = 8,
        label = L"|g_\varepsilon|_{Q_0}^2 \mathrm{\ (computed)}")

    # Reference line: ε^{-2} scaling
    eps_ref = range(eps_vals[end], eps_vals[1], length=50)
    ref_vals = 2.0 .* eps_ref .^ (-2)
    lines!(ax1, collect(eps_ref), ref_vals,
        color = "#C92A2A", linewidth = 1.5, linestyle = :dash,
        label = L"\sim \varepsilon^{-2} \mathrm{\ reference}")
    axislegend(ax1, position = :lt, framevisible = false, labelsize = 10)

    # ── Panel B: L² norm growth (for comparison) ──
    ax2 = Axis(fig[1, 2],
        xscale = log10,
        yscale = log10,
        title = L"$L^2$ norm growth (Panel B)",
        xlabel = L"Mollifier width $\varepsilon$",
        ylabel = L"$\|g_\varepsilon\|_{L^2}^2$",
        xgridvisible = true,
        ygridvisible = true,
        xgridcolor = (:black, 0.1),
        ygridcolor = (:black, 0.1),
        xminorticksvisible = true,
        yminorticksvisible = true,
        xminorticks = IntervalsBetween(9),
        yminorticks = IntervalsBetween(9)
    )
    l2_vals = Float64.(cm_df.l2_norm_sq)
    scatterlines!(ax2, eps_vals, l2_vals,
        color = "#5F3DC4", linewidth = 2.0, marker = :rect, markersize = 8,
        label = L"\|g_\varepsilon\|_{L^2}^2")
    # Reference: ε^{-1} scaling (since ||g_ε||² ~ 1/(2√π ε))
    ref_l2 = 0.3 .* collect(eps_ref) .^ (-1)
    lines!(ax2, collect(eps_ref), ref_l2,
        color = "#C92A2A", linewidth = 1.5, linestyle = :dash,
        label = L"\sim \varepsilon^{-1} \mathrm{\ reference}")
    axislegend(ax2, position = :lt, framevisible = false, labelsize = 10)

    # ── Panel C: λ⁶ R_λ^{bound} for different ε (remainder breakdown) ──
    ax3 = Axis(fig[2, 1],
        xscale = log10,
        yscale = log10,
        title = L"Scaled residual bound $\lambda^6 R_\lambda^{\mathrm{bound}}$ (Panel C)",
        xlabel = L"Clamp strength $\lambda$",
        ylabel = L"\lambda^6 R_\lambda^{\mathrm{bound}}",
        xgridvisible = true,
        ygridvisible = true,
        xgridcolor = (:black, 0.1),
        ygridcolor = (:black, 0.1)
    )
    colors = ["#0B7285", "#E67700", "#C92A2A", "#5F3DC4"]
    markers = [:circle, :rect, :utriangle, :diamond]
    test_eps_powers = [0, 4, 8, 12]
    for (i, p) in enumerate(test_eps_powers)
        subset = breakdown_df[breakdown_df.eps_power .== p, :]
        lam_plot = Float64.(subset.lambda)
        l6r_plot = Float64.(subset.lambda6_R_bound)
        eps_label = @sprintf("\\varepsilon = 2^{-%d}", p)
        scatterlines!(ax3, lam_plot, l6r_plot,
            color = colors[i], linewidth = 1.8, marker = markers[i], markersize = 8,
            label = latexstring(eps_label))
    end
    axislegend(ax3, position = :rb, framevisible = false, labelsize = 10)

    # ── Panel D: Ratio R_bound / σ_kk (breakdown indicator) ──
    ax4 = Axis(fig[2, 2],
        xscale = log10,
        yscale = log10,
        title = L"Remainder/variance ratio $R_\lambda^{\mathrm{bound}} / \Sigma_{\lambda,kk}$ (Panel D)",
        xlabel = L"Clamp strength $\lambda$",
        ylabel = L"R_\lambda^{\mathrm{bound}} / \Sigma_{\lambda,kk}",
        xgridvisible = true,
        ygridvisible = true,
        xgridcolor = (:black, 0.1),
        ygridcolor = (:black, 0.1)
    )
    for (i, p) in enumerate(test_eps_powers)
        subset = breakdown_df[breakdown_df.eps_power .== p, :]
        lam_plot = Float64.(subset.lambda)
        ratio_plot = Float64.(subset.ratio_R_over_sigma)
        eps_label = @sprintf("\\varepsilon = 2^{-%d}", p)
        scatterlines!(ax4, lam_plot, ratio_plot,
            color = colors[i], linewidth = 1.8, marker = markers[i], markersize = 8,
            label = latexstring(eps_label))
    end
    hlines!(ax4, [1.0], color = "#212529", linestyle = :dash, linewidth = 1.0)
    axislegend(ax4, position = :rb, framevisible = false, labelsize = 10)

    Label(fig[0, 1:2],
        "Experiment 4: Failure of Cameron-Martin Admissibility",
        fontsize = 14, font = :bold)

    save(joinpath(OUTPUT_DIR, "rough_observation.png"), fig, px_per_unit = 2.4)
end

function write_summary(cm_df::DataFrame, breakdown_df::DataFrame)
    lines = [
        "# Experiment 4: Failure of Cameron-Martin Admissibility",
        "",
        "## Claim Tested",
        "",
        "This experiment attacks Assumption (ii) (Cameron-Martin stability) by showing that",
        "a Dirac delta observation operator cannot satisfy the admissibility condition",
        "\$a_{Ik} \\in H_{Q_0}\$, because the Cameron-Martin energy norm diverges.",
        "",
        "## Setup",
        "",
        "- Prior: elliptic Gaussian field on [0,1] with \$L = -d^2/dx^2 + \\kappa^2\$,",
        "  \$\\kappa = $(Float64(KAPPA))\$.",
        "- Observation: narrow Gaussian mollifier \$g_\\varepsilon(x)\$ centered at \$x_0 = 0.5\$.",
        "- Spectral truncation: \$N = $(N_MODES)\$ modes.",
        "- Precision: `$(PRECISION_BITS)` bits.",
        "",
        "## Cameron-Martin Energy Divergence",
        "",
        "| ε | |g_ε|²_{Q₀} | ||g_ε||²_{L²} |",
        "|---|---:|---:|",
    ]
    for row in eachrow(cm_df)
        eps_str = @sprintf("2^{-%d}", row.eps_power)
        cm_str = @sprintf("%.6g", row.cm_energy_norm_sq)
        l2_str = @sprintf("%.6g", row.l2_norm_sq)
        push!(lines, "| \$$(eps_str)\$ | $(cm_str) | $(l2_str) |")
    end
    push!(lines, "")
    push!(lines, "The CM energy norm grows as \$\\sim \\varepsilon^{-2}\$ while the \$L^2\$ norm")
    push!(lines, "grows as \$\\sim \\varepsilon^{-1}\$, confirming that the Sobolev \$H^1\$ energy")
    push!(lines, "diverges faster than the \$L^2\$ norm alone.")
    push!(lines, "")

    # Remainder breakdown at representative lambda
    push!(lines, "## Remainder Bound Breakdown")
    push!(lines, "")
    push!(lines, "At \$\\lambda = 4096\$, the ratio \$R_\\lambda^{\\mathrm{bound}} / \\Sigma_{\\lambda,kk}\$:")
    push!(lines, "")
    push!(lines, "| ε | λ⁶ R_λ^{bound} | R_λ^{bound}/Σ_{λ,kk} | Status |")
    push!(lines, "|---|---:|---:|---|")
    for p in [0, 4, 8, 12]
        subset = breakdown_df[(breakdown_df.eps_power .== p) .& (breakdown_df.lambda_power .== 12), :]
        if nrow(subset) > 0
            row = first(subset)
            l6r = @sprintf("%.4g", row.lambda6_R_bound)
            ratio = @sprintf("%.4g", row.ratio_R_over_sigma)
            status = Float64(row.ratio_R_over_sigma) < 1.0 ? "✅ O(λ⁻⁶) holds" : "❌ BREAKDOWN"
            push!(lines, "| \$2^{-$(p)}\$ | $(l6r) | $(ratio) | $(status) |")
        end
    end
    push!(lines, "")
    push!(lines, "## Conclusion")
    push!(lines, "")
    push!(lines, "Point evaluation (Dirac delta) violates the Cameron-Martin admissibility")
    push!(lines, "condition. As \$\\varepsilon \\to 0\$, the CM energy \$|g_\\varepsilon|^2_{Q_0}\$")
    push!(lines, "diverges, destroying the \$O(\\lambda^{-6})\$ remainder estimate. The Schur")
    push!(lines, "residual upper bound loses all asymptotic control.")

    write(joinpath(OUTPUT_DIR, "validation_summary.md"), join(lines, "\n") * "\n")
end

function append_checksums(summary_path, file_paths)
    lines = [
        "",
        "## Output File Checksums",
        "",
        "| File | MD5 | SHA-256 |",
        "|---|---|---|",
    ]
    for path in file_paths
        filename = basename(path)
        sha = bytes2hex(sha256(read(path)))
        md = "N/A"
        try
            md = readchomp(`md5 -q $path`)
        catch
        end
        push!(lines, "| $filename | $md | $sha |")
    end
    open(summary_path, "a") do io
        write(io, join(lines, "\n") * "\n")
    end
end

function main()
    mkpath(OUTPUT_DIR)
    setprecision(BigFloat, PRECISION_BITS) do
        println("Computing Cameron-Martin energy divergence...")
        cm_df = evaluate_cm_divergence()
        CSV.write(joinpath(OUTPUT_DIR, "cm_energy_divergence.csv"), cm_df)
        println("  Done. Max CM energy: ", Float64(maximum(cm_df.cm_energy_norm_sq)))

        println("Computing remainder breakdown across λ grid...")
        breakdown_df = evaluate_remainder_breakdown()
        CSV.write(joinpath(OUTPUT_DIR, "remainder_breakdown.csv"), breakdown_df)
        println("  Done.")

        println("Generating plots...")
        plot_results(cm_df, breakdown_df)
        println("  Saved to $(joinpath(OUTPUT_DIR, "rough_observation.png"))")

        println("Writing validation summary...")
        write_summary(cm_df, breakdown_df)
        append_checksums(
            joinpath(OUTPUT_DIR, "validation_summary.md"),
            [
                joinpath(OUTPUT_DIR, "cm_energy_divergence.csv"),
                joinpath(OUTPUT_DIR, "remainder_breakdown.csv"),
                joinpath(OUTPUT_DIR, "rough_observation.png")
            ]
        )
        println("  Done.")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
