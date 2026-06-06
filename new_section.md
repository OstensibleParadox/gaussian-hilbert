由于 $C_\lambda \coloneqq \Sigma_{\lambda,II}$ 作为一个紧算子（且为迹类），它在无穷维空间上必然不存在有界全局逆。

通过 Tikhonov 正则化与 Galerkin 截断构造逼近序列，并利用 Anderson--Trapp 短路算子理论及自伴算子谱定理，可证：有限维或正则化计算给出的残差矩阵，在 Loewner 偏序（$\preceq$）下单调递增，并依算子范数收敛于残差 $R_{\lambda,E}$。这里的 Loewner 单调性只属于逼近序列。

\section{Finite-Rank Shorting and Exact Gaussian Finite-Part Cancellation}
\label{sec:finite_rank}

We have shown that singular Gaussian coordinate conditioning in a separable Hilbert space admits an inverse-free Schur normal form.  The main obstruction is the absence of a bounded inverse for compact covariance blocks.  The variational Cameron--Martin projection is the target-coordinate shorted-operator subtraction evaluated in the covariance-energy scale, and it yields the residual estimate needed for the hard-conditioning limit.  The rank-one construction reveals the basic inverse-free mechanism, but the finite-rank case requires an operator-valued residual and a matrix finite-part cancellation.  This section proves that finite-dimensional observed subspaces remain within the same Gaussian shorting theory.  The limitation is not finite rank, but infinite-rank observation and non-Gaussian conditioning.

\subsection{Matrix algebraic extension}
\label{sec:finite_rank:matrix_extension}

For an $m$-dimensional observed subspace $E=\operatorname{span}\{\phi_1,\ldots,\phi_m\}$, the scalar target block $\Sigma_{\lambda,kk}$ is replaced by the covariance matrix $\Sigma_{\lambda,EE}\in\mathbb{R}^{m\times m}$. The scalar calibration becomes a matrix Lyapunov balancing problem.

\begin{theorem}[Matrix Lyapunov Calibration]
  \label{thm:matrix_lyapunov}
  Let $E$ be an $m$-dimensional target subspace. Suppose the regularized drift and diffusion on $E$ are given by $A_{\lambda,E} = -\lambda \Lambda_E$ and $D_{\lambda,E} = \frac{1}{\lambda} \Delta_E$, where $\Lambda_E, \Delta_E \in \mathbb{R}^{m \times m}$ are strictly positive definite matrices. Let $\Sigma_{\lambda,EE}$ be the unique positive definite solution to the matrix Lyapunov equation:
  \begin{equation}
    \label{eq:discussion_matrix_lyapunov}
    A_{\lambda,E}\Sigma_{\lambda,EE}
    +\Sigma_{\lambda,EE}A_{\lambda,E}^{*}
    +D_{\lambda,E}=0.
  \end{equation}
  Then $\Sigma_{\lambda,EE}$ scales exactly as $\lambda^{-2}$. Specifically, the normalized matrix $\widetilde{\Sigma}_{EE} \coloneqq \lambda^2 \Sigma_{\lambda,EE}$ is independent of $\lambda$ and is the unique positive definite solution to the scale-free Lyapunov equation:
  \begin{equation}
    \label{eq:scale_free_lyapunov}
    \Lambda_E \widetilde{\Sigma}_{EE} + \widetilde{\Sigma}_{EE} \Lambda_E^* = \Delta_E.
  \end{equation}
\end{theorem}
\begin{proof}
  Substituting $A_{\lambda,E} = -\lambda \Lambda_E$ and $D_{\lambda,E} = \frac{1}{\lambda} \Delta_E$ into \eqref{eq:discussion_matrix_lyapunov} yields:
  \[
    -\lambda \Lambda_E \Sigma_{\lambda,EE} - \lambda \Sigma_{\lambda,EE} \Lambda_E^* + \frac{1}{\lambda} \Delta_E = 0.
  \]
  Multiplying the entire equation by $\lambda$ gives:
  \[
    \Lambda_E (\lambda^2 \Sigma_{\lambda,EE}) + (\lambda^2 \Sigma_{\lambda,EE}) \Lambda_E^* = \Delta_E.
  \]
  Defining $\widetilde{\Sigma}_{EE} \coloneqq \lambda^2 \Sigma_{\lambda,EE}$, we recover \eqref{eq:scale_free_lyapunov}. Since $\Lambda_E$ is positive definite, its spectrum lies in the open right half-plane, making $-\Lambda_E$ Hurwitz. Because $\Delta_E$ is positive definite, standard continuous Lyapunov theory guarantees that \eqref{eq:scale_free_lyapunov} admits a unique positive definite solution $\widetilde{\Sigma}_{EE}$. Therefore, the exact scaling $\Sigma_{\lambda,EE} = \lambda^{-2} \widetilde{\Sigma}_{EE}$ holds for all $\lambda > 0$, properly generalizing the scalar calibration condition $\Sigma_{\lambda,kk} = \varepsilon_0 / (2\lambda^2)$.
\end{proof}

\subsection{Operator-Valued Residuals and Loewner Approximation}
\label{sec:discussion:operator-valued-residuals}

For an $m$-dimensional target space $E$, the Schur residual must be a positive-semidefinite operator on $E$. Because the free-block covariance $C_\lambda \coloneqq \Sigma_{\lambda,II}$ is compact and lacks a bounded global inverse, the formal finite-dimensional matrix expression $R_{\lambda,E} = \Sigma_{\lambda,EI} C_\lambda^{-1} \Sigma_{\lambda,IE}$ is mathematically undefined. We must reconstruct the operator-valued residual directly through operator approximation sequences, tracking their convergence and monotonicity in the Loewner partial order. The Anderson--Trapp variational principle bridges the abstract range-inclusion theory to computationally constructive Galerkin and Tikhonov limits.

\begin{theorem}[Operator-Valued Schur Residual and Loewner Convergence]
  \label{thm:operator_residual}
  Let $\mathcal{H} = E \oplus \mathcal{H}_I$, with $E$ finite-dimensional. Let $\Sigma_\lambda$ be the non-negative trace-class covariance operator on $\mathcal{H}$ with blocks $\Sigma_{\lambda,EE}$, $\Sigma_{\lambda,EI}$, $\Sigma_{\lambda,IE}$, and $C_\lambda \coloneqq \Sigma_{\lambda,II}$. Assume the structural range inclusion $\operatorname{Range}(\Sigma_{\lambda,IE}) \subset \operatorname{Range}(C_\lambda^{1/2})$ and set
  \[
    H_{C_\lambda}:=\operatorname{Ran}(C_\lambda^{1/2}),
    \qquad
    \|h\|_{H_{C_\lambda}}
    :=
    \|C_\lambda^{-1/2}h\|_{\mathcal H_I}.
  \]
  Here $C_\lambda^{-1/2}$ denotes the RKHS inverse-image norm, not a bounded operator on the ambient Hilbert space. Define
  \[
    T_\lambda v := C_\lambda^{-1/2}\Sigma_{\lambda,IE}v,
    \qquad
    R_{\lambda,E}:=T_\lambda^*T_\lambda.
  \]
  Equivalently, $R_{\lambda,E}$ has quadratic forms
  \begin{equation}
    \label{eq:discussion_operator_residual_form}
    \langle v,R_{\lambda,E}v\rangle_E
    =
    \sup_{\langle u, C_\lambda u\rangle>0}
    \frac{|\langle u,\Sigma_{\lambda,IE}v\rangle|^2}
         {\langle u, C_\lambda u\rangle},
    \qquad v\in E.
  \end{equation}
  Then $R_{\lambda,E}$ is a well-defined bounded positive-semidefinite matrix on $E$. Furthermore, if $\|\Sigma_{\lambda,IE} v\|_{H_{C_\lambda}} \le M \lambda^{-3} \|v\|_E$ for all $v \in E$, then $\|R_{\lambda,E}\|_{\mathcal{L}(E)} = O(\lambda^{-6})$.

  Crucially, this exact residual is the norm-convergent limit of monotonically increasing sequences in the Loewner partial order:
  \begin{enumerate}
    \item \textbf{Galerkin Monotonicity:} Let $\Pi_n$ be a sequence of finite-rank orthogonal projections on $\mathcal{H}_I$ converging strongly to $I_{\mathcal{H}_I}$. The finite-dimensional Schur complements:
    \begin{equation}
      R_{\lambda,E}^{(n)} \coloneqq \Sigma_{\lambda,EI} \Pi_n (\Pi_n C_\lambda \Pi_n)^{\dagger} \Pi_n \Sigma_{\lambda,IE}
    \end{equation}
    form a Loewner-monotonic increasing sequence ($0 \le R_{\lambda,E}^{(n)} \preceq R_{\lambda,E}^{(n+1)} \preceq R_{\lambda,E}$), and converge to $R_{\lambda,E}$ in operator norm.
    \item \textbf{Tikhonov Monotonicity:} For any $\epsilon > 0$, the globally defined regularized residuals:
    \begin{equation}
      R_{\lambda,E}^{(\epsilon)} \coloneqq \Sigma_{\lambda,EI} (C_\lambda + \epsilon I)^{-1} \Sigma_{\lambda,IE}
    \end{equation}
    satisfy $R_{\lambda,E}^{(\epsilon_1)} \preceq R_{\lambda,E}^{(\epsilon_2)} \preceq R_{\lambda,E}$ whenever $\epsilon_1 > \epsilon_2 > 0$, and converge to $R_{\lambda,E}$ in operator norm as $\epsilon \searrow 0$.
  \end{enumerate}
\end{theorem}
\begin{proof}
  By the assumed range inclusion, for any $v \in E$, the vector $y_v \coloneqq \Sigma_{\lambda,IE} v$ lies in $H_{C_\lambda} = \operatorname{Ran}(C_\lambda^{1/2})$. Let $w_v=C_\lambda^{-1/2}y_v$ denote its minimal inverse image in the RKHS norm.

  The quadratic form \eqref{eq:discussion_operator_residual_form} evaluates exactly to $\|w_v\|_{\mathcal H_I}^2=\|y_v\|_{H_{C_\lambda}}^2$. This implies $\langle v, R_{\lambda,E} v \rangle_E = \|T_\lambda v\|_{\mathcal H_I}^2$, making $R_{\lambda,E}$ a unique positive-semidefinite matrix on $E$. The energy bound guarantees $\|R_{\lambda,E}\|_{\mathcal{L}(E)} = O(\lambda^{-6})$.

  For the Galerkin monotonicity, evaluating the finite-dimensional matrix inversion yields exactly the restricted supremum:
  \[
    \langle v, R_{\lambda,E}^{(n)} v \rangle_E = \sup_{u \in \operatorname{Range}(\Pi_n), \langle u, C_\lambda u \rangle > 0} \frac{|\langle u, y_v \rangle|^2}{\langle u, C_\lambda u \rangle}.
  \]
  Because $\operatorname{Range}(\Pi_n) \subset \operatorname{Range}(\Pi_{n+1})$, the supremum is taken over expanding sets, directly implying the Loewner monotonicity $\langle v, R_{\lambda,E}^{(n)} v \rangle_E \le \langle v, R_{\lambda,E}^{(n+1)} v \rangle_E \le \|w_v\|^2$. Since $E$ is finite-dimensional, any Loewner-monotonic bounded sequence of symmetric matrices converges in operator norm. The strong convergence of $\Pi_n \to I_{\mathcal{H}_I}$ guarantees the limit is the full supremum $R_{\lambda,E}$.

  For the Tikhonov monotonicity, note that $R_{\lambda,E}^{(\epsilon)}$ can be rewritten using $w_v$:
  \[
    \langle v, R_{\lambda,E}^{(\epsilon)} v \rangle_E = \langle y_v, (C_\lambda + \epsilon I)^{-1} y_v \rangle = \langle w_v, C_\lambda (C_\lambda + \epsilon I)^{-1} w_v \rangle.
  \]
  By the spectral theorem applied to the compact positive operator $C_\lambda$, the operator family $f_\epsilon(C_\lambda) \coloneqq C_\lambda(C_\lambda + \epsilon I)^{-1}$ is monotonically increasing as $\epsilon \searrow 0$ and converges strongly to the orthogonal projection onto $\overline{\operatorname{Range}(C_\lambda)}$. Since $w_v \in \overline{\operatorname{Range}(C_\lambda)}$, the quadratic form converges monotonically to $\|w_v\|^2$, establishing both the Loewner order $R_{\lambda,E}^{(\epsilon_1)} \preceq R_{\lambda,E}^{(\epsilon_2)}$ and the norm convergence to $R_{\lambda,E}$.
\end{proof}

\subsection{General observation operators and misalignment}
\label{sec:discussion:misalignment}

The scalar proof in this paper uses a block decomposition aligned with a single observed coordinate. However, typical finite-rank observation operators $L:\mathcal H\to\mathbb R^m$ select a subspace $E$ that does not commute with the free prior covariance, meaning the target-free split is not an eigenbasis split. This geometric misalignment generates cross-coupling that requires a precise resolvent analysis.

\begin{theorem}[Misaligned Covariance Coupling]
  \label{thm:misaligned_coupling}
  Let $P_E$ be an orthogonal projection onto an $m$-dimensional observed subspace $E$, and let $\mathcal{H}_I = E^\perp$. Let the regularized dynamics be decomposed such that the target block matches \Cref{thm:matrix_lyapunov} with $A_{\lambda,E} = -\lambda \Lambda_E$, the cross-noise is zero ($D_{\lambda,IE} = 0$), and the off-diagonal drift coupling is denoted by $A_{\lambda,IE}$.
  Assume that, for some $\lambda_0>0$,
  \[
    \sup_{\lambda\ge\lambda_0}
    \|A_{\lambda,IE}\|_{\mathcal L(E,H_{Q_0})}<\infty.
  \]
  Assume also that the restricted drift $A_{II}$ generates an exponentially stable semigroup on $H_{Q_0}$ and that there exist $K,\gamma>0$ such that
  \[
    \|e^{-t\Lambda_E^*}\|_{\mathcal L(E)}
    \le
    K e^{-\gamma t},
    \qquad t\ge0.
  \]
  Then the cross-covariance satisfies the exact identity:
  \begin{equation}
    \label{eq:misaligned_resolvent_identity}
    \Sigma_{\lambda,IE} = \int_0^\infty e^{t A_{II}} A_{\lambda,IE} \Sigma_{\lambda,EE} e^{t A_{\lambda,E}^*} dt.
  \end{equation}
  Furthermore, the misaligned cross-covariance energy scales as $\|\Sigma_{\lambda,IE} v\|_{H_{Q_0}} = O(\lambda^{-3})$ for any $v \in E$.
\end{theorem}
\begin{proof}
  The off-diagonal block of the global Lyapunov equation $\mathcal{A}^{(\lambda)} \Sigma_\lambda + \Sigma_\lambda (\mathcal{A}^{(\lambda)})^* + \mathcal{D}^{(\lambda)} = 0$ yields the Sylvester equation:
  \[
    A_{II} \Sigma_{\lambda,IE} + \Sigma_{\lambda,IE} A_{\lambda,E}^* + A_{\lambda,IE} \Sigma_{\lambda,EE} = 0.
  \]
  Since $A_{II}$ generates an exponentially stable semigroup $e^{t A_{II}}$ and the assumed bound controls $e^{tA_{\lambda,E}^*}=e^{-\lambda t\Lambda_E^*}$, the unique solution is given by the integral representation \eqref{eq:misaligned_resolvent_identity}.
  
  To bound the Cameron--Martin norm, fix $v \in E$. We have $e^{t A_{\lambda,E}^*} v = e^{-\lambda t \Lambda_E^*} v$. By \Cref{thm:matrix_lyapunov}, $\Sigma_{\lambda,EE} = O(\lambda^{-2})$. By the uniform coupling assumption, choose $C$ such that $\|A_{\lambda,IE} w\|_{H_{Q_0}} \le C \|w\|_E$ for all $w \in E$ and all $\lambda\ge\lambda_0$. Using the uniform exponential stability of $e^{t A_{II}}$ on $H_{Q_0}$, namely $\|e^{t A_{II}}\|_{\mathcal{L}(H_{Q_0})} \le M e^{-\omega t}$ for some $M \ge 1, \omega > 0$, we estimate:
  \begin{align*}
    \|\Sigma_{\lambda,IE} v\|_{H_{Q_0}} 
    &\le \int_0^\infty \|e^{t A_{II}}\|_{\mathcal{L}(H_{Q_0})} \|A_{\lambda,IE} \Sigma_{\lambda,EE} e^{-\lambda t \Lambda_E^*} v\|_{H_{Q_0}} dt \\
    &\le \int_0^\infty M e^{-\omega t} C \|\Sigma_{\lambda,EE}\|_{\mathcal{L}(E)} \|e^{-\lambda t \Lambda_E^*}\|_{\mathcal{L}(E)} \|v\|_E dt.
  \end{align*}
  Using $\|e^{-\lambda t\Lambda_E^*}\|_{\mathcal{L}(E)} \le K e^{-\lambda\gamma t}$, we obtain:
  \begin{align*}
    \|\Sigma_{\lambda,IE} v\|_{H_{Q_0}} 
    &\le M C K \|\Sigma_{\lambda,EE}\|_{\mathcal{L}(E)} \|v\|_E \int_0^\infty e^{-(\omega + \lambda \gamma) t} dt \\
    &= M C K \|\Sigma_{\lambda,EE}\|_{\mathcal{L}(E)} \|v\|_E \frac{1}{\omega + \lambda \gamma}.
  \end{align*}
  Since $\|\Sigma_{\lambda,EE}\|_{\mathcal{L}(E)} = O(\lambda^{-2})$ and $(\omega + \lambda \gamma)^{-1} = O(\lambda^{-1})$, we conclude that $\|\Sigma_{\lambda,IE} v\|_{H_{Q_0}} = O(\lambda^{-3})$. Combined with \Cref{thm:operator_residual}, this verifies that the shorted-operator residual maintains the $O(\lambda^{-6})$ scaling even in the presence of observation misalignment.
\end{proof}

\subsection{Multidimensional Gaussian Renormalized Finite-Part Functional}
\label{sec:discussion:multidim_functional}

With the finite-rank covariance bounds in \Cref{thm:matrix_lyapunov,thm:operator_residual,thm:misaligned_coupling} in place, the rank-one finite-part functional from \Cref{sec:gaussian_zero} becomes the prototype for the matrix case. In the multidimensional case, the finite-part functional must be evaluated in the metric of the residual matrix, and the KL divergence subtraction requires operator trace identities rather than scalar logarithm expansions. We now re-derive the exact cancellation using the Weinstein--Aronszajn identity.

The scalar renormalized functional can be fully generalized to an $m$-dimensional observation subspace $E$. Let the exact evidence be $X_E = x_E^* \in \mathbb{R}^m$. The true Schur residual matrix is $S_{\lambda,E} \coloneqq \Sigma_{\lambda,EE} - R_{\lambda,E} \in \mathbb{R}^{m \times m}$, where $R_{\lambda,E}$ is the shorted-operator residual from \Cref{thm:operator_residual}. Since $R_{\lambda,E}^{(n)}\nearrow R_{\lambda,E}$ in the Loewner order, the Galerkin residual matrices
\[
S_{\lambda,E}^{(n)}
=
\Sigma_{\lambda,EE}-R_{\lambda,E}^{(n)}
\]
decrease in the Loewner order to
\[
S_{\lambda,E}
=
\Sigma_{\lambda,EE}-R_{\lambda,E}.
\]
They are not claimed to decrease strictly.  Positivity of the limiting
matrix is guaranteed by the condition
\[
\rho_{\lambda,E}
=
\Sigma_{\lambda,EE}^{-1/2}
R_{\lambda,E}
\Sigma_{\lambda,EE}^{-1/2}
<I.
\]

\begin{definition}[Multidimensional Finite-Part Functional]
  \label{def:multidim_unified_action}
  Whenever $S_{\lambda,E} \succ 0$ and $\mu_{\text{obs},I}$ is equivalent to $\mu_{\lambda,I}$, the multidimensional Gaussian renormalized finite-part functional is defined by evaluating the conditional residual in the metric of the Schur LMMSE matrix:
  \begin{equation}
    \label{eq:multidim_unified_action}
    \begin{aligned}
      \mathfrak{J}_{\lambda,E}(\mu_{\mathrm{obs}}; \mu_\lambda) \coloneqq &\KL(\mu_{\mathrm{obs},I} \parallel \mu_{\lambda,I}) \\
      &+ \frac{1}{2} \E_{\mu_{\mathrm{obs}}}\left[ \left( X_E - \mu_{E|I}(X_I) \right)^T S_{\lambda,E}^{-1} \left( X_E - \mu_{E|I}(X_I) \right) \right],
    \end{aligned}
  \end{equation}
  where $\mu_{E|I}(X_I) \coloneqq m_{\lambda,E} + W_{\Sigma_{\lambda,IE}}(X_I - m_{\lambda,I})$ is the measurable LMMSE predictor on $\mathbb{R}^m$, with the Wiener integral evaluated componentwise.
\end{definition}

To establish the multidimensional limit, we define the dimensionless coupling matrix $\rho_{\lambda,E} \coloneqq \Sigma_{\lambda,EE}^{-1/2} R_{\lambda,E} \Sigma_{\lambda,EE}^{-1/2} \in \mathbb{R}^{m \times m}$. By the energy bound in \Cref{thm:operator_residual}, $\|\rho_{\lambda,E}\| = O(\lambda^{-4})$.

\begin{theorem}[Exact Matrix Trace and Determinant Cancellation]
  \label{thm:multidim_cancellation}
  For any evidence level $t \in \mathbb{R}^m$, let $\mu_{\lambda}^t$ be the weakly continuous conditional law on $X_E = t$. The infinite-dimensional KL divergence and the matrix-valued residual expectation undergo exact cancellation of all cross-coupling traces and quadratic forms, yielding:
  \begin{equation}
    \label{eq:multidim_exact_cancellation}
    \mathfrak{J}_{\lambda,E}(\mu_\lambda^t; \mu_\lambda) = \frac{1}{2} (t - m_{\lambda,E})^T \Sigma_{\lambda,EE}^{-1} (t - m_{\lambda,E}) - \frac{1}{2}\log\det(I-\rho_{\lambda,E}).
  \end{equation}
\end{theorem}
\begin{proof}
  Let $C_\lambda \coloneqq \Sigma_{\lambda,II}$. Under $\mu_\lambda^t$, the free coordinates follow $\mathcal{N}(\mu_{I|E}(t), \Sigma_{I|E})$. By the Feldman--H\'{a}jek formula for trace-class perturbations, evaluating the KL divergence relies on the operator $K \coloneqq C_\lambda^{-1/2} \Sigma_{\lambda,IE} \Sigma_{\lambda,EE}^{-1/2}$. By the Weinstein--Aronszajn identity, the non-zero spectrum of the infinite-dimensional operator $KK^*$ matches that of the $m \times m$ matrix $K^*K = \rho_{\lambda,E}$. This reduces the infinite-dimensional Fredholm determinant and trace exactly to the matrix manifold:
  \[
    \KL(\mu_{\lambda,I}^t \parallel \mu_{\lambda,I}) = \frac{1}{2} \left[ (t - m_{\lambda,E})^T \Sigma_{\lambda,EE}^{-1} R_{\lambda,E} \Sigma_{\lambda,EE}^{-1} (t - m_{\lambda,E}) - \log\det(I - \rho_{\lambda,E}) - \tr(\rho_{\lambda,E}) \right].
  \]
  For the residual expectation term, let $Z \coloneqq t - \mu_{E|I}(X_I)$. Using block matrix identities, the conditional mean is $\E_{\mu_\lambda^t}[Z] = S_{\lambda,E} \Sigma_{\lambda,EE}^{-1} (t - m_{\lambda,E})$ and the conditional variance is $\Var_{\mu_\lambda^t}(Z) = R_{\lambda,E} \Sigma_{\lambda,EE}^{-1} S_{\lambda,E}$. Expanding the quadratic form $\frac{1}{2}\E[Z^T S_{\lambda,E}^{-1} Z]$ gives:
  \[
    \frac{1}{2} \left[ (t - m_{\lambda,E})^T \Sigma_{\lambda,EE}^{-1} (t - m_{\lambda,E}) - (t - m_{\lambda,E})^T \Sigma_{\lambda,EE}^{-1} R_{\lambda,E} \Sigma_{\lambda,EE}^{-1} (t - m_{\lambda,E}) + \tr(\rho_{\lambda,E}) \right].
  \]
  Summing the KL divergence and the quadratic residual form exactly cancels both the $\pm \frac{1}{2}\tr(\rho_{\lambda,E})$ trace anomalies and the cross-covariance residual matrices, isolating \eqref{eq:multidim_exact_cancellation}.
\end{proof}

\begin{theorem}[Multidimensional Gaussian Finite-Part Regularization]
  \label{thm:multidim_regularization}
  Assume the $m$-dimensional regularized drift on $E$ is centered at the observed value:
  \[
    A_{\lambda,E}(m_{\lambda,E}-x_E^*)+b_E=0,
    \qquad
    A_{\lambda,E}=-\lambda\Lambda_E.
  \]
  Defining the normalized displacement $\delta_E \coloneqq \Lambda_E^{-1} b_E \in \mathbb{R}^m$, the mean deterministic offset satisfies
  \[
    m_{\lambda,E}-x_E^*
    =
    \lambda^{-1}\Lambda_E^{-1}b_E
    =
    \lambda^{-1}\delta_E.
  \]
  
  As $\lambda \to \infty$, the multidimensional renormalized finite-part functional converges to a well-defined bounded limit governed by the scale-free Lyapunov solution $\widetilde{\Sigma}_{EE}$:
  \begin{equation}
    \lim_{\lambda\to\infty} \mathfrak{J}_{\lambda,E}(\mu_{\mathrm{obs}}; \mu_\lambda) = \frac{1}{2} \delta_E^T \widetilde{\Sigma}_{EE}^{-1} \delta_E.
  \end{equation}
  The approximation residual satisfies $\mathfrak{J}_{\lambda,E} = \frac{1}{2} \delta_E^T \widetilde{\Sigma}_{EE}^{-1} \delta_E + O(\lambda^{-4})$.
\end{theorem}
\begin{proof}
  Substitute $t = x_E^*$ and $t - m_{\lambda,E} = -\lambda^{-1} \delta_E$ into \Cref{thm:multidim_cancellation}. Under the canonical scaling,
  \[
    \Sigma_{\lambda,EE}
    =
    \lambda^{-2}\widetilde\Sigma_{EE},
    \qquad
    t-m_{\lambda,E}
    =
    -\lambda^{-1}\delta_E,
    \qquad
    \|\rho_{\lambda,E}\|=O(\lambda^{-4}).
  \]
  Therefore the quadratic term yields exactly
  \[
    \frac12
    \delta_E^T
    \widetilde\Sigma_{EE}^{-1}
    \delta_E,
  \]
  and the determinant term satisfies $\log\det(I-\rho_{\lambda,E})=O(\lambda^{-4})$. Hence
  \[
    \mathfrak{J}_{\lambda,E}
    =
    \frac12
    \delta_E^T
    \widetilde\Sigma_{EE}^{-1}
    \delta_E
    +
    O(\lambda^{-4}).
  \]
\end{proof}

\subsection{Non-Gaussian and dynamical extensions}
\label{sec:discussion:non-gaussian}

The result is therefore deliberately limited to covariance and Gaussian-conditioning questions in the finite-rank aligned settings.  It does not assert a general non-Gaussian conditioning theory.  For non-Gaussian stationary processes, the projection residual $S_{\lambda,E}=\Sigma_{\lambda,EE}-R_{\lambda,E}$ still has a second-moment interpretation, but Radon disintegration, finite-part functionals, and measure-class behavior require separate hypotheses \cite{bogachev2007measure,kallenberg2021foundations}.  Dynamical relaxation after conditioning, or transport distances between post-conditioning laws, belongs to a different problem.  The present contribution is the operator-theoretic construction of the singular Gaussian conditioning residual itself.
