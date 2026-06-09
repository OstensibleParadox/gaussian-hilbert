"""
Experiment 5: Failure of Dual Scaling / Constant Diffusion.

This stress test attacks the canonical scaling condition (Proposition 2,
Theorem 3 in the paper) by fixing D_{λ,E} = D₀ (constant background noise)
instead of the canonical D_{λ,E} = Δ_E/λ.

Under the canonical scaling:
  - A_{λ,E} = -λΛ_E,  D_{λ,E} = Δ_E/λ
  - Lyapunov: Σ_{λ,EE} = λ⁻² Σ̃_{EE}         (exact λ⁻² scaling)

Under constant diffusion (this experiment):
  - A_{λ,E} = -λΛ_E,  D_{λ,E} = D₀ = const
  - Lyapunov: Σ_{λ,EE} = D₀/(2λΛ_E) ~ λ⁻¹   (wrong scaling!)

Consequences of the broken scaling:
  1. Σ_{EE} ~ λ⁻¹ instead of λ⁻²
  2. Cross-covariance energy: Σ_{IE} ~ λ⁻² instead of λ⁻³
  3. Shorted residual: R_{λ,E} ~ λ⁻⁴ instead of λ⁻⁶
  4. Fredholm cancellation fails: logdet and trace terms have
     mismatched asymptotic orders

The experiment uses a 2D target subspace to show the matrix-level
failure of the Fredholm trace-anomaly cancellation.
"""

using CairoMakie
using LaTeXStrings
using DataFrames
using CSV
using Printf
using SHA
using LinearAlgebra

const OUTPUT_DIR = joinpath(@__DIR__, "results")
const PRECISION_BITS = 256
const LAMBDA_POWERS = 0:12

# ── Model Parameters ──
# Drift on target E (2D): A_{λ,E} = -λ Λ_E
# Λ_E = diag(1, 2)  (two target modes with spectral gaps 1 and 2)
# Canonical diffusion: Δ_E = diag(1, 1), so D_{λ,E} = Δ_E/λ
# Constant diffusion: D₀ = diag(1, 1)

# Free block (1D for simplicity): a scalar background mode
# A_II = -ω (scalar, ω = spectral gap of background)
# D_II = d_I (scalar, background noise)
# Off-diagonal coupling: A_{IE} = [c₁, c₂]' (background driven by target)

const LAMBDA_E = [1.0, 2.0]   # target drift eigenvalues
const DELTA_E  = [1.0, 1.0]   # canonical diffusion eigenvalues
const D0_E     = [1.0, 1.0]   # constant diffusion (same matrix, but not scaled)

const OMEGA_BG  = 3.0          # background drift eigenvalue
const D_BG      = 1.0          # background noise
const COUPLING  = [0.5, 0.3]   # A_{IE} coupling vector (2D → 1D)

"""
    lyapunov_target(lam, lambda_e, d_e)

Solve the decoupled target Lyapunov equation for a diagonal system:
    A_{λ,E} Σ_{EE} + Σ_{EE} A_{λ,E}' + D_{λ,E} = 0
where A_{λ,E} = -λ diag(Λ_E) and D_{λ,E} = diag(d_e).

Returns the diagonal target covariance Σ_{EE} = diag(d_e ./ (2λΛ_E)).
"""
function lyapunov_target(lam::BigFloat, lambda_e::Vector{Float64}, d_e::Vector{BigFloat})
    m = length(lambda_e)
    sigma_ee = zeros(BigFloat, m, m)
    for i in 1:m
        sigma_ee[i, i] = d_e[i] / (BigFloat(2) * lam * BigFloat(lambda_e[i]))
    end
    return sigma_ee
end

"""
    cross_covariance(lam, sigma_ee, lambda_e, coupling, omega)

Compute the cross-covariance Σ_{IE} from the off-diagonal Sylvester equation:
    A_II Σ_{IE} + Σ_{IE} A_{λ,E}' + A_{IE} Σ_{EE} = 0

For a 1D background with A_II = -ω (scalar) and A_{λ,E} = -λ diag(Λ_E):
    Σ_{IE}[1,j] = coupling[j] * Σ_{EE}[j,j] / (ω + λ Λ_E[j])
"""
function cross_covariance(lam::BigFloat, sigma_ee::Matrix{BigFloat},
                          lambda_e::Vector{Float64}, coupling::Vector{Float64},
                          omega::Float64)
    m = length(lambda_e)
    sigma_ie = zeros(BigFloat, 1, m)  # 1×m since background is 1D
    for j in 1:m
        sigma_ie[1, j] = BigFloat(coupling[j]) * sigma_ee[j, j] /
                         (BigFloat(omega) + lam * BigFloat(lambda_e[j]))
    end
    return sigma_ie
end

"""
    background_covariance(lam, sigma_ee, sigma_ie, lambda_e, coupling, omega, d_bg)

Compute the background covariance Σ_{II} from the (1,1) block of the
Lyapunov equation. For 1D background:
    -2ω Σ_II + 2 coupling' Σ_{IE}' + d_bg = 0
    Σ_II = (d_bg + 2 Σ_{IE} coupling) / (2ω)
(Here coupling' * Σ_{IE}' gives the forcing from target fluctuations.)
"""
function background_covariance(lam::BigFloat, sigma_ee::Matrix{BigFloat},
                                sigma_ie::Matrix{BigFloat},
                                coupling::Vector{Float64}, omega::Float64,
                                d_bg::Float64)
    forcing = BigFloat(d_bg)
    for j in eachindex(coupling)
        forcing += BigFloat(2) * BigFloat(coupling[j]) * sigma_ie[1, j]
    end
    return forcing / (BigFloat(2) * BigFloat(omega))
end

"""
    schur_residual(sigma_ee, sigma_ie, sigma_ii)

Compute the Schur residual matrix:
    R_{λ,E} = Σ_{IE}' Σ_{II}^{-1} Σ_{IE}

In operator-theoretic terms this is the finite-dimensional evaluation of
the Anderson-Trapp supremum (exact when the free block is 1D).
"""
function schur_residual(sigma_ee::Matrix{BigFloat}, sigma_ie::Matrix{BigFloat},
                        sigma_ii::BigFloat)
    m = size(sigma_ee, 1)
    R = zeros(BigFloat, m, m)
    for i in 1:m, j in 1:m
        R[i, j] = sigma_ie[1, i] * sigma_ie[1, j] / sigma_ii
    end
    return R
end

"""
    evaluate_scaling(lam, diffusion_family)

Evaluate all covariance blocks and diagnostic quantities for a given λ
and diffusion family ("canonical" or "constant").

Returns a named tuple with all diagnostic values.
"""
function evaluate_scaling(lam::BigFloat, family::String)
    m = length(LAMBDA_E)

    # Target diffusion
    if family == "canonical"
        d_e = [BigFloat(DELTA_E[i]) / lam for i in 1:m]
    elseif family == "constant"
        d_e = [BigFloat(D0_E[i]) for i in 1:m]
    else
        error("Unknown family: $family")
    end

    # Solve Lyapunov blocks
    sigma_ee = lyapunov_target(lam, LAMBDA_E, d_e)
    sigma_ie = cross_covariance(lam, sigma_ee, LAMBDA_E, COUPLING, OMEGA_BG)
    sigma_ii = background_covariance(lam, sigma_ee, sigma_ie, COUPLING, OMEGA_BG, D_BG)

    # Schur residual
    R = schur_residual(sigma_ee, sigma_ie, sigma_ii)
    S = sigma_ee - R  # Shorted Schur complement

    # Diagnostic quantities
    trace_sigma_ee = tr(sigma_ee)
    trace_R = tr(R)
    trace_S = tr(S)

    # Normalized coupling ratio: ρ = Σ_{EE}^{-1/2} R Σ_{EE}^{-1/2}
    # For diagonal Σ_{EE}: ρ[i,j] = R[i,j] / sqrt(Σ_{EE}[i,i] * Σ_{EE}[j,j])
    # ||ρ|| = tr(ρ) for rank-1 PSD (opnorm = trace = unique eigenvalue)
    # Computed entirely in BigFloat to avoid Float64 precision loss.
    rho = zeros(BigFloat, m, m)
    for i in 1:m, j in 1:m
        rho[i, j] = R[i, j] / sqrt(sigma_ee[i, i] * sigma_ee[j, j])
    end
    norm_rho = tr(rho)  # BigFloat trace; equals opnorm for rank-1 PSD

    # Fredholm functional components (for the trace anomaly)
    # logdet(I - ρ) and tr(ρ) should cancel at the same order under canonical scaling
    # Compute in BigFloat: -log(1 - ρ_jj) for each diagonal entry
    trace_rho = norm_rho  # same as tr(ρ)
    logdet_term = zero(BigFloat)
    for j in 1:m
        logdet_term -= log(one(BigFloat) - rho[j, j])
    end
    # Trace anomaly: logdet(I - ρ) + tr(ρ) ≈ 0.5 tr(ρ²) + O(ρ³)
    trace_anomaly = logdet_term - trace_rho

    return (
        family = family,
        lambda = lam,
        sigma_ee_11 = sigma_ee[1, 1],
        sigma_ee_22 = sigma_ee[2, 2],
        trace_sigma_ee = trace_sigma_ee,
        lam_sq_trace_sigma = lam^2 * trace_sigma_ee,
        trace_R = trace_R,
        lam6_trace_R = lam^6 * trace_R,
        trace_S = trace_S,
        norm_rho = norm_rho,
        logdet_term = logdet_term,
        trace_rho = trace_rho,
        trace_anomaly = trace_anomaly,
        lam4_trace_anomaly = lam^4 * trace_anomaly,
    )
end

"""
    asymptotic_rho_probe(lam)

Compute ||ρ|| = tr(ρ) and its asymptotic decomposition at a given λ, all in BigFloat.

Exact formula (rank-1 background, diagonal Σ_{EE}):
    tr(ρ) = Σⱼ Σ_{IE,j}² / (Σ_{II} Σ_{EE,jj})
          = Σⱼ cⱼ² Σ_{EE,jj} / ((ω + λΛⱼ)² Σ_{II})

Substituting Σ_{EE,jj} = Δⱼ/(2λ²Λⱼ) for canonical scaling:
    tr(ρ) = Σⱼ cⱼ² Δⱼ / (2λ² Λⱼ (ω + λΛⱼ)² Σ_{II})

Leading-order asymptotic (λ → ∞):
    tr(ρ)_∞ = (ω / (d_bg λ⁴)) Σⱼ cⱼ² Δⱼ / Λⱼ³   =  O(λ⁻⁴)

The exact/asymptotic ratio reveals the sub-leading ω/(λΛⱼ) correction.
"""
function asymptotic_rho_probe(lam::BigFloat)
    m = length(LAMBDA_E)
    omega = BigFloat(OMEGA_BG)
    d_bg = BigFloat(D_BG)

    # Canonical diffusion
    d_e = [BigFloat(DELTA_E[i]) / lam for i in 1:m]
    sigma_ee = lyapunov_target(lam, LAMBDA_E, d_e)
    sigma_ie = cross_covariance(lam, sigma_ee, LAMBDA_E, COUPLING, OMEGA_BG)
    sigma_ii = background_covariance(lam, sigma_ee, sigma_ie, COUPLING, OMEGA_BG, D_BG)

    # Exact ||ρ|| = tr(ρ) (rank-1 PSD ⟹ opnorm = trace = unique eigenvalue)
    norm_rho_exact = zero(BigFloat)
    for j in 1:m
        norm_rho_exact += sigma_ie[1, j]^2 / (sigma_ii * sigma_ee[j, j])
    end

    # Leading-order asymptotic at λ → ∞:
    #   Σ_{EE,jj} → Δⱼ/(2λ²Λⱼ), Σ_{II} → d_bg/(2ω), (ω+λΛⱼ) → λΛⱼ
    #   tr(ρ)_∞ = (2ω/d_bg) Σⱼ cⱼ²Δⱼ / (2λ⁴Λⱼ³) = (ω/(d_bg λ⁴)) Σⱼ cⱼ²Δⱼ/Λⱼ³
    sigma_ii_inf = d_bg / (BigFloat(2) * omega)
    norm_rho_leading = zero(BigFloat)
    for j in 1:m
        c_j = BigFloat(COUPLING[j])
        L_j = BigFloat(LAMBDA_E[j])
        D_j = BigFloat(DELTA_E[j])
        # Each term: cⱼ² · Σ_{EE,jj,∞} / ((λΛⱼ)² · σ_{II,∞})
        #          = cⱼ² · Dⱼ/(2λ²Λⱼ) / (λ²Λⱼ² · d_bg/(2ω))
        #          = cⱼ² Dⱼ ω / (d_bg λ⁴ Λⱼ³)
        norm_rho_leading += c_j^2 * D_j * omega / (d_bg * lam^4 * L_j^3)
    end

    # Next-order correction from (1 + ω/(λΛⱼ))⁻² ≈ 1 - 2ω/(λΛⱼ) + ...
    # Weighted correction coefficient for the O(1/λ) relative term:
    #   tr(ρ)/tr(ρ)_∞ ≈ 1 + c₁/λ + O(1/λ²)
    #   c₁ = -2ω · (Σ cⱼ²Δⱼ/Λⱼ⁴) / (Σ cⱼ²Δⱼ/Λⱼ³)
    numer_baseline = zero(BigFloat)
    numer_correction = zero(BigFloat)
    for j in 1:m
        c_j = BigFloat(COUPLING[j])
        L_j = BigFloat(LAMBDA_E[j])
        D_j = BigFloat(DELTA_E[j])
        numer_baseline += c_j^2 * D_j / L_j^3
        numer_correction += c_j^2 * D_j / L_j^4
    end
    correction_coeff = -BigFloat(2) * omega * numer_correction / numer_baseline

    return (
        lambda = lam,
        norm_rho_exact = norm_rho_exact,
        norm_rho_leading = norm_rho_leading,
        ratio_exact_over_leading = norm_rho_exact / norm_rho_leading,
        sigma_ii_exact = sigma_ii,
        sigma_ii_inf = sigma_ii_inf,
        sigma_ii_ratio = sigma_ii / sigma_ii_inf,
        lam4_rho = lam^4 * norm_rho_exact,
        correction_coeff = correction_coeff,  # multiply by 1/λ for the O(1/λ) correction
        predicted_ratio = one(BigFloat) + correction_coeff / lam,
    )
end


function plot_results(df::DataFrame)
    fig = Figure(size=(900, 720))

    canonical = df[df.family .== "canonical", :]
    constant = df[df.family .== "constant", :]

    # ── Panel A: λ² tr(Σ_{EE}) scaling ──
    ax1 = Axis(fig[1, 1],
        xscale = log10,
        yscale = log10,
        title = L"Target covariance scaling $\lambda^2 \mathrm{tr}(\Sigma_{EE})$ (Panel A)",
        xlabel = L"Clamp strength $\lambda$",
        ylabel = L"\lambda^2 \mathrm{tr}(\Sigma_{\lambda,EE})",
        xgridvisible = true,
        ygridvisible = true,
        xgridcolor = (:black, 0.1),
        ygridcolor = (:black, 0.1)
    )
    # Canonical: λ² tr(Σ_EE) → const
    lam_c = Float64.(canonical.lambda)
    scatterlines!(ax1, lam_c, Float64.(canonical.lam_sq_trace_sigma),
        color = "#0B7285", linewidth = 2.0, marker = :circle, markersize = 8,
        label = L"Canonical $D_\lambda = \Delta_E / \lambda$")
    # Constant: λ² tr(Σ_EE) ~ λ (diverges)
    lam_k = Float64.(constant.lambda)
    scatterlines!(ax1, lam_k, Float64.(constant.lam_sq_trace_sigma),
        color = "#C92A2A", linewidth = 2.0, marker = :rect, markersize = 8,
        label = L"Constant $D_\lambda = D_0$")
    hlines!(ax1, [0.75], color = "#212529", linestyle = :dash, linewidth = 1.0)
    axislegend(ax1, position = :lt, framevisible = false, labelsize = 10)

    # ── Panel B: Normalized coupling ||ρ|| ──
    ax2 = Axis(fig[1, 2],
        xscale = log10,
        yscale = log10,
        title = L"Normalized coupling $\|\rho_{\lambda,E}\|$ (Panel B)",
        xlabel = L"Clamp strength $\lambda$",
        ylabel = L"\|\rho_{\lambda,E}\|",
        xgridvisible = true,
        ygridvisible = true,
        xgridcolor = (:black, 0.1),
        ygridcolor = (:black, 0.1)
    )
    scatterlines!(ax2, lam_c, Float64.(canonical.norm_rho),
        color = "#0B7285", linewidth = 2.0, marker = :circle, markersize = 8,
        label = L"Canonical: $\|\rho\| \sim \lambda^{-4}$")
    scatterlines!(ax2, lam_k, Float64.(constant.norm_rho),
        color = "#C92A2A", linewidth = 2.0, marker = :rect, markersize = 8,
        label = L"Constant: $\|\rho\| \sim \lambda^{-2}$")
    axislegend(ax2, position = :rt, framevisible = false, labelsize = 10)

    # ── Panel C: λ⁶ tr(R) ──
    ax3 = Axis(fig[2, 1],
        xscale = log10,
        yscale = log10,
        title = L"Scaled residual $\lambda^6 \mathrm{tr}(R_{\lambda,E})$ (Panel C)",
        xlabel = L"Clamp strength $\lambda$",
        ylabel = L"\lambda^6 \mathrm{tr}(R_{\lambda,E})",
        xgridvisible = true,
        ygridvisible = true,
        xgridcolor = (:black, 0.1),
        ygridcolor = (:black, 0.1)
    )
    scatterlines!(ax3, lam_c, Float64.(canonical.lam6_trace_R),
        color = "#0B7285", linewidth = 2.0, marker = :circle, markersize = 8,
        label = L"Canonical: $\lambda^6 R \to \mathrm{const}$")
    # For constant diffusion, λ⁶ R ~ λ² (diverges), so we plot it
    scatterlines!(ax3, lam_k, Float64.(constant.lam6_trace_R),
        color = "#C92A2A", linewidth = 2.0, marker = :rect, markersize = 8,
        label = L"Constant: $\lambda^6 R \to \infty$")
    axislegend(ax3, position = :lt, framevisible = false, labelsize = 10)

    # ── Panel D: Trace anomaly ──
    ax4 = Axis(fig[2, 2],
        xscale = log10,
        yscale = log10,
        title = L"Fredholm trace anomaly $\lambda^4 \cdot \mathrm{anomaly}$ (Panel D)",
        xlabel = L"Clamp strength $\lambda$",
        ylabel = L"\lambda^4 |\mathrm{logdet}(I-\rho) - \mathrm{tr}(\rho)|",
        xgridvisible = true,
        ygridvisible = true,
        xgridcolor = (:black, 0.1),
        ygridcolor = (:black, 0.1)
    )
    # For canonical: λ⁴ anomaly → const (the cancellation works)
    canonical_anomaly = [Float64(abs(x)) for x in canonical.lam4_trace_anomaly]
    canonical_anomaly = [max(x, 1e-20) for x in canonical_anomaly]
    scatterlines!(ax4, lam_c, canonical_anomaly,
        color = "#0B7285", linewidth = 2.0, marker = :circle, markersize = 8,
        label = L"Canonical: $\lambda^4 \cdot \mathrm{anomaly} \to \mathrm{const}$")
    # For constant: λ⁴ anomaly diverges (wrong order)
    constant_anomaly = [Float64(abs(x)) for x in constant.lam4_trace_anomaly]
    constant_anomaly = [max(x, 1e-20) for x in constant_anomaly]
    scatterlines!(ax4, lam_k, constant_anomaly,
        color = "#C92A2A", linewidth = 2.0, marker = :rect, markersize = 8,
        label = L"Constant: diverges")
    axislegend(ax4, position = :lt, framevisible = false, labelsize = 10)

    Label(fig[0, 1:2],
        "Experiment 5: Failure of Dual Scaling (Constant Diffusion)",
        fontsize = 14, font = :bold)

    save(joinpath(OUTPUT_DIR, "constant_diffusion.png"), fig, px_per_unit = 2.4)
end

function write_summary(df::DataFrame)
    canonical = df[df.family .== "canonical", :]
    constant = df[df.family .== "constant", :]
    can_final = last(sort(canonical, :lambda))
    con_final = last(sort(constant, :lambda))

    lines = [
        "# Experiment 5: Failure of Dual Scaling (Constant Diffusion)",
        "",
        "## Claim Tested",
        "",
        "The canonical diffusion scaling \$D_{\\lambda,E} = \\Delta_E/\\lambda\$ is",
        "**structurally necessary**, not merely a cosmetic assumption. Replacing it",
        "with a constant diffusion \$D_{\\lambda,E} = D_0\$ breaks the exact",
        "\$\\lambda^{-2}\$ covariance scaling, which in turn destroys the Fredholm",
        "trace-anomaly cancellation.",
        "",
        "## Model Setup",
        "",
        "- 2D target subspace \$E\$ with drift eigenvalues \$\\Lambda_E = \\mathrm{diag}($(LAMBDA_E[1]), $(LAMBDA_E[2]))\$.",
        "- 1D scalar background mode with spectral gap \$\\omega = $(OMEGA_BG)\$.",
        "- Off-diagonal coupling \$A_{IE} = [$(COUPLING[1]), $(COUPLING[2])]'\$.",
        "- Precision: `$(PRECISION_BITS)` bits.",
        "",
        "## Scaling Comparison at λ = $(Float64(can_final.lambda))",
        "",
        "| Quantity | Canonical | Constant | Ratio |",
        "|---|---:|---:|---:|",
        @sprintf("| λ² tr(Σ_{EE}) | %.8f | %.4g | — |", can_final.lam_sq_trace_sigma, con_final.lam_sq_trace_sigma),
        @sprintf("| λ⁶ tr(R) | %.8f | %.4g | — |", can_final.lam6_trace_R, con_final.lam6_trace_R),
        @sprintf("| ||ρ|| | %.4g | %.4g | — |", can_final.norm_rho, con_final.norm_rho),
        @sprintf("| λ⁴ anomaly | %.4g | %.4g | — |", abs(can_final.lam4_trace_anomaly), abs(con_final.lam4_trace_anomaly)),
        "",
        "## Diagnosis",
        "",
        "Under constant diffusion:",
        "- Σ_{EE} ~ λ⁻¹ (instead of λ⁻²): the target variance decays too slowly.",
        "- Cross-covariance energy ~ λ⁻² (instead of λ⁻³): the coupling is too strong.",
        "- Shorted residual R ~ λ⁻⁴ (instead of λ⁻⁶): remainder is too large.",
        "- Normalized coupling ||ρ|| ~ λ⁻² (instead of λ⁻⁴): the trace anomaly terms",
        "  (logdet and trace) diverge relative to the intended O(λ⁻⁴) cancellation.",
        "",
        "The Fredholm cancellation fundamentally fails because the cross-term",
        "\$\\Sigma_{IE}\\Sigma_{EE}^{-1}\\Sigma_{EI}\$ and the logdet/trace terms",
        "no longer share the same asymptotic order.",
    ]
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
        records = NamedTuple[]

        println("Evaluating canonical scaling family...")
        for p in LAMBDA_POWERS
            lam = BigFloat(2)^p
            push!(records, evaluate_scaling(lam, "canonical"))
        end
        println("  Done.")

        println("Evaluating constant diffusion family...")
        for p in LAMBDA_POWERS
            lam = BigFloat(2)^p
            push!(records, evaluate_scaling(lam, "constant"))
        end
        println("  Done.")

        df = DataFrame(records)
        CSV.write(joinpath(OUTPUT_DIR, "constant_diffusion_data.csv"), df)

        println("Generating plots...")
        plot_results(df)
        println("  Saved to $(joinpath(OUTPUT_DIR, "constant_diffusion.png"))")

        println("Writing validation summary...")
        write_summary(df)
        append_checksums(
            joinpath(OUTPUT_DIR, "validation_summary.md"),
            [
                joinpath(OUTPUT_DIR, "constant_diffusion_data.csv"),
                joinpath(OUTPUT_DIR, "constant_diffusion.png")
            ]
        )
        println("  Done.")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
