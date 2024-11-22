# Matlab source

This folder contains both the implementations of the various walk-based Laplacians and the codes to generate the numerical experiments contained in the paper.

## Walk-Laplacian Implementations

| code |   | Implementation  |
|------|---|-----------------|
|`kwalklaplacian.m`| $\mathbb{L}(f) = \text{diag}(f(A)\mathbf{1}) - f(A)$ | Implementation of the walk-based Laplacians |
|`buildkatzbasedlaplacian.m` | $\mathbb{L}^{NBT}(\alpha) = (1-\alpha^2)(\text{diag}( (I - \alpha A - \alpha^2(I-D))^{-1}\mathbf{1}) - (I - \alpha A - \alpha^2(I-D))^{-1} )$ | Implementation of the NBT Katz based Laplacian, the matrix-vector product routine is implemented in `katzmatvec.m` |

