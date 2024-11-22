# Matlab source

This folder contains both the implementations of the various walk-based Laplacians and the codes to generate the numerical experiments contained in the paper.

## Walk-Laplacian Implementations

<table>
<tr>
<th>Code</th>
<th></th>
<th>Implementation</th>
</tr>
<tr>
<td> kwalklaplacian.m </td>
<td>
$\mathbb{L}(f) = \text{diag}(f(A)\mathbf{1}) - f(A)$
</td>
<td>
Implementation of the walk-based Laplacians
</td>
</tr>
<tr>
<td>
buildkatzbasedlaplacian.m
</td>
<td>
$\mathbb{L}^{NBT}(\alpha) = (1-\alpha^2)(\text{diag}( (I - \alpha A - \alpha^2(I-D))^{-1}\mathbf{1}) - (I - \alpha A - \alpha^2(I-D))^{-1} )$
</td>
<td>
Implementation of the NBT Katz based Laplacian, the matrix-vector product routine is implemented in katzmatvec.m.
</td>
</tr>
<tr>
<td>
buildphibasedlaplacian.m
</td>
<td> $\mathbb{L}^{\text{NBT}}(\mathbf{c})$
</td>
<td>
Implementation of the NBT $\phi$-based Laplacian, the matrix-vector product routine si implemented in phimatvec.m. The code uses Algorithm 919.
</td>
</tr>
</table>



