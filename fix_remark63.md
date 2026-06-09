# Remark 6.3 修正笔记：暂态谱失衡与一致良态保护

## 1. 原始假说：暂态谱失衡（Transient Spectral Imbalance）

$$R_E^{(k)} = \operatorname{diag}(1-\epsilon, \ldots, 1-\epsilon, 0, \ldots, 0) \implies \kappa(S^{(k)}) = 1/\epsilon$$

这是一个合法的线性代数事实：各向异性 Galerkin 解析 + 高相关（$\rho$ 接近 $I$）= 暂态条件数尖峰。

## 2. 因果谬误的彻底厘清

* ❌ 伪因果：Loewner 非严格单调 $\implies$ 暂态条件数尖峰。
* ✅ 真因果：极高相关性（$\rho \approx I$） + 各向异性解析 $\implies$ 条件数尖峰。

非严格单调只是因为 Galerkin 空间扩充时"没有包含交叉协方差的全部满秩信息"，这是一个纯粹的几何/运动学现象（Kinematic constraint）。而条件数爆炸需要巨大的能量抽离，这是一个动力学现象（Dynamics）。既然系统天然是极弱耦合（$\rho = O(\lambda^{-4})$），真因果链的动力学齿轮就被彻底拆掉了。

## 3. Galerkin 逼近的真实理论生态

$S_{\lambda,E}^{(n)}$ 的 Galerkin 收敛是用来**证明极限存在**，不是作为计算 procedure。

引入 Anderson–Trapp 变分和 Galerkin 嵌套子空间，纯粹是为了在严谨的算子代数意义下，利用单调序列定理（Monotone Sequence Theorem）构造性地证明极限残差在 Hilbert 空间上是良定的。它根本不需要"严格单调"，只要"非递减且有上界"就大功告成了。

## 4. 暂态谱失衡的一般理论地位

在一般的多重网格法（Multigrid）、区域分解（Domain Decomposition）或强相关贝叶斯推断中，如果目标变量几乎可以被背景变量完全解释（即高度相关，$\rho \approx I$），那么随着残差被逐步消解，必然有 $\epsilon \to 0$。

* 此时叠加各向异性解析（即单调性不严格，某些方向更新了，某些停滞），中间态矩阵 $S^{(k)}$ 的谱就会被拉扯到极致，条件数狂飙到 $1/\epsilon$。

## 5. Singular Conditioning 框架下的无情镇压

但这在 Singular Conditioning 框架下，被底层的渐近物理量无情镇压了。

因为强迫正则化（Regularization clamp $\lambda$）切断目标与背景纠缠的速度，远快于压制目标自身方差的速度：

* 目标块方差基数：$\Sigma_{\lambda,EE} = O(\lambda^{-2})$
* 交叉协方差能量：$R_{\lambda,E} = O(\lambda^{-6})$
* 无量纲耦合度：$\|\rho_{\lambda,E}\| = O(\lambda^{-4}) \to 0$

这意味着，对于任意的中间态 $S_{\lambda,E}^{(n)} = \Sigma_{\lambda,EE} - R_{\lambda,E}^{(n)}$，由于 $0 \le R_{\lambda,E}^{(n)} \le R_{\lambda,E}$，对其做白化归一化后必然有：

$$I \ge \Sigma_{\lambda,EE}^{-1/2} S_{\lambda,E}^{(n)} \Sigma_{\lambda,EE}^{-1/2} \ge I - \rho_{\lambda,E}$$

因为 $\|\rho_{\lambda,E}\| = O(\lambda^{-4})$，中间矩阵的所有特征值被死死地夹在 $[1 - O(\lambda^{-4}), 1]$ 之间。

正如上界：

$$\kappa(S_{\lambda,E}^{(n)}) \le \kappa(\Sigma_{\lambda,EE}) \cdot \frac{1}{1 - O(\lambda^{-4})}$$

尖峰连"凸起"的机会都没有。系统不仅在极限处是良态的，在整条 Galerkin 逼近路径上，都受到了极度严密的**一致良态**（Uniform well-conditioning）保护。
