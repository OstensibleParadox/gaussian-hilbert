# Experiment 5: Failure of Dual Scaling (Constant Diffusion)

## Claim Tested

The canonical diffusion scaling $D_{\lambda,E} = \Delta_E/\lambda$ is
**structurally necessary**, not merely a cosmetic assumption. Replacing it
with a constant diffusion $D_{\lambda,E} = D_0$ breaks the exact
$\lambda^{-2}$ covariance scaling, which in turn destroys the Fredholm
trace-anomaly cancellation.

## Model Setup

- 2D target subspace $E$ with drift eigenvalues $\Lambda_E = \mathrm{diag}(1.0, 2.0)$.
- 1D scalar background mode with spectral gap $\omega = 3.0$.
- Off-diagonal coupling $A_{IE} = [0.5, 0.3]'$.
- Precision: `256` bits.

## Scaling Comparison at λ = 4096.0

| Quantity | Canonical | Constant | Ratio |
|---|---:|---:|---:|
| λ² tr(Σ_{EE}) | 0.75000000 | 3072 | — |
| λ⁶ tr(R) | 0.38288261 | 6.424e+06 | — |
| ||ρ|| | 2.78e-15 | 1.139e-11 | — |
| λ⁴ anomaly | 9.983e-16 | 1.675e-08 | — |

## Diagnosis

Under constant diffusion:
- Σ_{EE} ~ λ⁻¹ (instead of λ⁻²): the target variance decays too slowly.
- Cross-covariance energy ~ λ⁻² (instead of λ⁻³): the coupling is too strong.
- Shorted residual R ~ λ⁻⁴ (instead of λ⁻⁶): remainder is too large.
- Normalized coupling ||ρ|| ~ λ⁻² (instead of λ⁻⁴): the trace anomaly terms
  (logdet and trace) diverge relative to the intended O(λ⁻⁴) cancellation.

The Fredholm cancellation fundamentally fails because the cross-term
$\Sigma_{IE}\Sigma_{EE}^{-1}\Sigma_{EI}$ and the logdet/trace terms
no longer share the same asymptotic order.

## Output File Checksums

| File | MD5 | SHA-256 |
|---|---|---|
| constant_diffusion_data.csv | 87041024f439138e158eac09961f7021 | fe84ba43ec893289921bf696d42f38376efb28c9b226958b4e271d7ac0d06693 |
| constant_diffusion.png | 348435d7d670715d8a6b625489a4d3f1 | 6f6bced0e9ed019d913f976102fc528dcaf6735983f22e079ad6f66ffcf7870b |
