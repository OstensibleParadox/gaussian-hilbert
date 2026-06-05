# Singular Gaussian Conditioning in Hilbert Space: Schur Residuals, Fisher Scaling, and Cameron--Martin Norms

This repository contains the LaTeX source code and numerical experiments for the paper:
**"Singular Gaussian Conditioning in Hilbert Space: Schur Residuals, Fisher Scaling, and Cameron--Martin Norms"**

Author: **Lucia Yizi Zhang** (Independent Researcher)
Email: [lucia.zhang.research@outlook.com](mailto:lucia.zhang.research@outlook.com)

---

## Repository Structure

The numerical validation of the paper's main theorems is organized into experiments corresponding to sections in **Appendix D** of the paper:

```
paper1_pr_it/
├── README.md                      # This documentation file
├── experiments/                   # Numerical experiments source code
│   ├── exp1_truncation_convergence/ # Appendix D.1: Truncation Convergence
│   │   ├── validate_truncation_convergence.jl
│   │   ├── Project.toml           # Julia project dependencies
│   │   ├── Manifest.toml          # Julia environment state
│   │   └── results/               # CSV outputs and summary logs
│   ├── exp2_covariance_calibration/ # Appendix D.2: Exact Covariance Calibration Ablation
│   │   ├── gaussian_clamped_action_ablation.jl
│   │   ├── Project.toml
│   │   ├── Manifest.toml
│   │   └── results/
│   └── exp3_three_node_chain/     # Appendix D.3: Three-Node Gaussian Chain
│       ├── three_node_gaussian_chain.jl
│       ├── Project.toml
│       ├── Manifest.toml
│       └── results/
└── paper/                         # Publication source files
    ├── paper1_merged.tex          # Main LaTeX source document
    ├── paper1_merged_blind.tex    # Anonymous submission LaTeX source document
    ├── appendix_merged.tex        # Unified appendix LaTeX source document
    ├── references.bib             # Bibliography database
    ├── paper1_merged.pdf          # Compiled publication PDF
    ├── paper1_merged_blind.pdf    # Compiled anonymous publication PDF
    └── figures/                   # Compiled figures for the paper
        ├── truncation_convergence_validation.png
        ├── ablation_combined.png
        ├── three_node_gaussian_chain.png
        ├── dichotomy_geometries.png
        └── mahalanobis_well_convergence.png
```

---

## Experiment Details & Claims Tested

### 1. Appendix D.1: Truncation Convergence (`experiments/exp1_truncation_convergence/`)
* **Claim:** Evaluates the Galerkin truncation convergence behavior predicted by Theorem 1 against a high-resolution reference covariance on an admissible infinite-dimensional Ornstein–Uhlenbeck witness.
* **Output:** Generates `truncation_convergence.csv`, `lambda_schur_scaling.csv`, the plot `truncation_convergence_validation.png` (copied to the paper's figure folder), and `validation_summary.md`.

### 2. Appendix D.2: Exact Covariance Calibration Ablation (`experiments/exp2_covariance_calibration/`)
* **Claim:** Validates the exact covariance calibration property in its decoupled reference channel. Tests four target-diffusion families (Canonical, Slow, Fast, Perturbed) using 256-bit `BigFloat` arithmetic, showing that only the canonical family satisfies the exact finite-$\lambda$ condition.
* **Output:** Generates `ablation_results.csv`, `ablation_summary.md`, and the plot `ablation_combined.png`.

### 3. Appendix D.3: Three-Node Gaussian Chain (`experiments/exp3_three_node_chain/`)
* **Claim:** Evaluates the Gaussian renormalized finite-part functional $\mathfrak{J}_\lambda$ on a three-node linear SDE ($X_1 \to X_2 \to X_3$) across three regimes: decoupled-centered (exact zero), coupled-centered ($O(\lambda^{-4})$ decay), and coupled-off-centered (converges to $1$).
* **Output:** Generates `three_node_results.csv`, `three_node_summary.md`, and the plot `three_node_gaussian_chain.png`.

### 4. Additional Appendix D Diagnostics
The unified manuscript also incorporates:
* **Appendix D.4: Numerical Precision Audit:** Compares IEEE 754 `Float64` and 256-bit Julia `BigFloat` arithmetic, showing that block-decoupling suppresses roundoff propagation down to $10^{-17}$ relative error.
* **Appendix D.5: Dichotomy of Dual Geometries:** Validates the simultaneous boundary survival of the 0th-order flat metric and the $O(\lambda^{-2})$ collapse of the unrescaled horizontal metric (Figure: `dichotomy_geometries.png`).
* **Appendix D.6: Orthogonal Independence:** Numerically verifies that the 0th-order normal geometry (unperturbed noise scale) and 1st-order horizontal lift (off-diagonal Cameron--Martin norm) vary independently.
* **Appendix D.7: Convergence to the Mahalanobis Normal Well:** Tracks the $O(\lambda^{-4})$ convergence of the blow-up geometry contours to the limiting Mahalanobis well (Figure: `mahalanobis_well_convergence.png`).

---

## Cryptographic Verification & Reproducibility

To guarantee reproducibility and prevent data corruption, every experiment automatically computes and appends the **MD5** and **SHA-256** checksums of its output CSV data and generated figures.

These checksums are logged in the `results/validation_summary.md` file inside each experiment directory. After reproducing the experiments locally, you can calculate the checksums of your local files and verify them against the committed records in this repository to confirm bit-perfect replication.

---

## Reproduction Instructions

### Prerequisites
* **Julia:** Version 1.9 or higher is recommended.
* **Libraries:** The environment is pre-configured with `Project.toml` and `Manifest.toml` files.

### Step-by-Step Execution

First, navigate to the paper repository:
```bash
cd paper1_pr_it
```

#### Run Experiment 1 (Truncation Convergence)
```bash
# Instantiate environment and run script
julia --project=experiments/exp1_truncation_convergence -e 'using Pkg; Pkg.instantiate()'
julia --project=experiments/exp1_truncation_convergence experiments/exp1_truncation_convergence/validate_truncation_convergence.jl
```

#### Run Experiment 6 (Exact Calibration Ablation)
```bash
julia --project=experiments/exp2_covariance_calibration -e 'using Pkg; Pkg.instantiate()'
julia --project=experiments/exp2_covariance_calibration experiments/exp2_covariance_calibration/gaussian_clamped_action_ablation.jl
```

#### Run Experiment 7 (Three-Node Gaussian Chain)
```bash
julia --project=experiments/exp3_three_node_chain -e 'using Pkg; Pkg.instantiate()'
julia --project=experiments/exp3_three_node_chain experiments/exp3_three_node_chain/three_node_gaussian_chain.jl
```

All generated tables and plot assets will be written to the respective `results/` subdirectory in each experiment folder. The plot images are automatically generated at double resolution (`px_per_unit = 2`) for inclusion in the publication PDF.
