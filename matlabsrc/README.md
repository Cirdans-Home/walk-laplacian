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

## Randomized Average Return Probability Estimator

To study the behavior of the diffusion process induced by different Laplacian operators we look at the average return probability

```math
\hat{p}_0(t) = \frac{1}{N} \sum_{i=1}^{N} [P(t)]_{i,i} = \frac{1}{N} \text{tr}\left( \exp(- t L ) \right) = \frac{1}{N} \sum_{i=1}^{N} \exp(- \lambda_i( L ) t),
```

The implementation of a stochastic estimator for its approximation is given in the function `xnystraceexp.m`:

```matlab
function [T,AVT,ERRT] = xnystraceexp(A,m,type,t,nt,varargin)
%XNYSTRACE-EXP randomized average return probability estimate
%   A adjacency matrix of the graph
%   m number of matrix vector product to be used in the estimate
%   type type of dynamic to be studied, admissible types are:
%       - 'LAPLACIAN' Standard Laplacian dynamic y'(t) = L^T y(t)
%   t = [t0,t1] time interval for the estimate, if longer uses the value
%   nt number of time points in the [t0,t1] interval, if t is longer than 2
%       takes the length of t as value
end
```

The implementation is based on the XTrace estimator, if you use this routine also mention:
> Epperly, Ethan N., Joel A. Tropp, and Robert J. Webber. "Xtrace: Making the most of every sample in stochastic trace estimation." SIAM Journal on Matrix Analysis and Applications 45.1 (2024): 1-23.

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
<tr>
<td> generate_tree.py </td>
<td>

```python
# This function generates a tree with a given number of nodes and a given branching factor using
# the NetworkX library. The script has to be executed through the MATLAB engine API in order to
# be able to call the function from MATLAB and return the adjacency matrix of the generated tree
# as a MATLAB sparse matrix.
```

</td>
</tr>
</table>

## Numerical Examples

These are the scripts needed to replicate the numerical experiments contained in the paper.

1. Graph with a trap `trap_example.m`. This example depens on the 919 algorithm, see the `get919.sh` script in the parent folder.
2. Random-tree example `tree_example.m`. This example depens on the 919 algorithm, see the `get919.sh` script in the parent folder.
   To run it you need to have a working interface between Python and MATLAB. On a Linux based system this can be achieved by running
   on the bash:
   ```bash
   python3 -m venv <path>
   source <path>/bin/activate
   pip install networkx numpy scipy
   ```
   Then going back to MATLAB and running
   ```matlab
   pyenv('Version','<path>/bin/python3','ExecutionMode','OutOfProcess')
   ```
   
> [!TIP]
> You can change the environment only before executing any Python script. If MATLAB complains, you should just restart it and execute
> the previous command as first thing.
3. Average Return Probabilities experiments
   1. The test in `averagereturnprobabilities.m` runs the example with the full computation, hence it is extremely taxing on both memory
      and time. In our case we run it on a cluster node under SLURM. The bash script `run_averagereturn.sh` has an example with our cluster
      configuratio. It should be easily adaptable to other system. The experiment can then be run by executing the command
      ```bash
      sbatch run_averagereturn.sh
      ```
      outputs are then saved as `.mat` files that can be used to produce the figures employing the `batchaveragereturn_plot.m` file.
   2. Examples with the stochastic estimator are contained in the files
      - `testnystrom.m` which contains the small scale example on the Karate complex network.
      - `batchnystrom.m` that runs the example on the larger networks. As for the other example on the larger network, we have run it
      in a node of the cluster, and an example of the batch script is contained in `runbatchnystrom.sh` that can be executed via
      ```bash
      sbatch runbatchnystrom.sh
      ```
      
> [!WARNING]
> The figures for the paper have been generated by means of the [matlab2tikz](https://github.com/matlab2tikz/matlab2tikz) package.
> It is not necessary to have it installed to run tests, the codes will complain about its existence saying they cannot produce figures, 
> but the experiments will still run.
