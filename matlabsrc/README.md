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

## Path-Laplacians

In order to perform some of the comparisons we need an implementation of the path-Laplacian (actually the transformed Path-Laplacian). This is implemented in the
file `pathlap.m`.

```matlab
function Lpath = pathlap(G,type,parameter)
    %PATHLAP Builds the transformed k-Path Laplacian of the given type
    %   INPUT: G graph object
    %          type = 'exp' or 'power'
    %          parameter = scaling parameter
    %   OUTPUT: Lpath transformed path Laplacian
end
```

If you use it, please cite the original works where this has been introduced:

> Estrada, Ernesto, et al. "Path Laplacian operators and superdiffusive processes on graphs. I. One-dimensional case." Linear Algebra and its Applications 523 (2017): 307-334.

> Estrada, Ernesto, et al. "Path Laplacian operators and superdiffusive processes on graphs. II. Two-dimensional lattice." Linear Algebra and its Applications 555 (2018): 373-397.


## Auxiliary routines

The code requires some auxiliary routines for parameter determinations, e.g., the $\alpha$ for the Katz-based Laplacian of NBT type, or the construction of Krylov spaces. These are briefly reported here.

<table>
<tr>
<th>
Code
</th>
<th>
Description
</th>
</tr>
<tr>
<td>producealphas.m</td>
<td>

```matlab
function [alphavec,rhoz] = producealphas(A,k,check)
%PRODUCEALPHAS Prints the value of alpha to be used in the Fortran test
%for evaluating the linear system solution routine.
```

See the listalphas.m file for its usage.
</td>
</tr>
</table>
