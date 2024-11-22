# Matlab source

This folder contains both the implementations of the various walk-based Laplacians and the codes to generate the numerical experiments contained in the paper.

## Walk-Laplacian Implementations

| code |   | Implementation  |
|------|---|-----------------|
|`kwalklaplacian.m`| $\mathbb{L}(f) = \text{diag}(f(A)\mathbf{1}) - f(A)$ | Implementation of the walk-based Laplacians |
|`buildkatzbasedlaplacian.m` | $\mathbb{L}^{NBT}(\alpha) = (1-\alpha^2)(\text{diag}( (I - \alpha A - \alpha^2(I-D))^{-1}\mathbf{1}) - (I - \alpha A - \alpha^2(I-D))^{-1} )$ | Implementation of the NBT Katz based Laplacian, the matrix-vector product routine is implemented in `katzmatvec.m` |
|`buildphibasedlaplacian.m`| $\mathbb{L}^{\text{NBT}}(\mathbf{c}) = $ | \text{diag}\left( \begin{bmatrix} O & I \end{bmatrix} \varphi_1(Z) \begin{bmatrix} \mathbf{d} \\ A^2\mathbf{1} - \mathbf{d} \end{bmatrix}  \right) -  \begin{bmatrix} O & I \end{bmatrix} \varphi_1(Z) \begin{bmatrix} A \\ A^2 - D \end{bmatrix}, | Implementation of the NBT $\varphi$-based Laplacian, the matrix-vector product routine si implemented in `phimatvec.m`. The code uses Algorithm 919.|

