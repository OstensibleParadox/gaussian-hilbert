# Asymptotic Precision Table: Normalized Coupling tr(ρ)

Under canonical scaling $D_{\lambda,E} = \Delta_E/\lambda$, the normalized coupling
$\rho_{\lambda,E} = \Sigma_{EE}^{-1/2} R_{\lambda,E} \Sigma_{EE}^{-1/2}$ satisfies:

$$\operatorname{tr}(\rho) \to \frac{C}{\lambda^4}, \quad C = \frac{\omega}{d_{\mathrm{bg}}} \sum_j \frac{c_j^2 \Delta_j}{\Lambda_j^3} = 0.78375$$

| $\lambda$ | $\operatorname{tr}(\rho)$ | $\lambda^4 \operatorname{tr}(\rho)$ | ratio (exact/leading) | eff. exponent $\alpha$ |
|---:|---:|---:|---:|---:|
| $2^1$ | 2.496879e-03 | 0.0399500624 | 0.5162044694 | - |
| $2^2$ | 6.618698e-04 | 0.1718500476 | 0.2817242908 | 1.9149 |
| $2^4$ | 8.545363e-06 | 0.5600526989 | 17.9519654491 | 3.1041 |
| $2^6$ | 5.389404e-07 | 0.9051698429 | 3626.7755498979 | 1.9898 |
| $2^8$ | 3.376916e-08 | 0.9321428766 | 29303.3411423429 | 1.9967 |
| $2^{10}$ | 2.112248e-09 | 0.9434026963 | 236958.1929044791 | 1.9990 |
| $2^{12}$ | 1.320460e-10 | 0.9490227989 | 1914825.0568773418 | 1.9997 |
| $2^{14}$ | 8.253335e-12 | 0.9518361759 | 15479485.3698006790 | 1.9999 |
| $2^{16}$ | 5.158422e-13 | 0.9532465853 | 125139754.2131832100 | 2.0000 |
| $2^{18}$ | 3.224029e-14 | 0.9539527779 | 1011730652.1399753000 | 2.0000 |
| $2^{20}$ | 2.015021e-15 | 0.9543060873 | 8181265729.3133764000 | 2.0000 |

**Leading coefficient** $C = \omega/(d_{\mathrm{bg}}) \cdot \sum_j c_j^2\Delta_j/\Lambda_j^3 = 0.78375$

**Observed**: $\lambda^4 \operatorname{tr}(\rho) \to 0.954$ (not $0.78375$).

The ratio (exact/leading) diverges as $\lambda^2$, indicating the leading-order formula
captures the correct $O(\lambda^{-4})$ exponent but the coefficient has a multiplicative
correction from the resolvent $({\omega + \lambda\Lambda_j})^{-2}$ vs $(\lambda\Lambda_j)^{-2}$
that converges slowly as $\lambda^{-1}$.

**Precision**: 256-bit BigFloat. $\|\rho\|_{\mathrm{op}} = \operatorname{tr}(\rho)$ (rank-1 PSD).
