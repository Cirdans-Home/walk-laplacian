# Walk based Laplacians for Modeling Diffusion on Complex Networks

This repository contains the accompanying code for the article "Walk based Laplacians for Modeling Diffusion on Complex Networks".

You can **clone the main repository** and initialize all submodules at once with this command:
```bash
git clone --recurse-submodules https://github.com/Cirdans-Home/walk-laplacian.git
```
the `--recurse-submodules` tells Git to automatically initialize and update each submodule. 

>[!TIP]
> If you’ve already cloned the repository without initializing submodules, you can initialize and update them afterward:
> ```bash
> git submodule update --init --recursive
> ```

## Collaborators

- Francesca Arrigo [:email:](mailto:francesca.arrigo@strath.ac.uk)
- Fabio Durastante [:email:](mailto:fabio.durastante@unipi.it) [:link:](https://fdurastante.github.io/)

## External codes

The code in this repository makes use of a number of external codes. In particular, it contains as Git submodules the [:link: PSCToolkit library](https://psctoolkit.github.io/) which is needed to repeat some of the experiments involving the solution of linear systems related to the Katz centrality-based Laplacian. Information about configuration and installation can be found in the READMEs contained in the :open_file_folder: configure and :open_file_folder: install folders.

The other needed codes are the Rational-Krylov Toolbox [RKToolbox](http://guettel.com/rktoolbox/) that can be obtained and installed by executing in the MATLAB shell the command:
```matlab
unzip('http://guettel.com/rktoolbox/rktoolbox.zip'); 
cd('rktoolbox'); addpath(fullfile(cd)); savepath
```
and the [:link: Algorithm 919](https://dl.acm.org/doi/abs/10.1145/2168773.2168781#sec-supp).
The `get919.sh` scripts downloads it and put it into the `matlabsrc` folder. If you are not running
under Linux, you can download the zip file from the supplementary materials of the paper at the previous link.

> [!IMPORTANT] 
> If you use the code contained in this repository and, in particular, the features inherited from external codes, please cite also the reference works for them.
> ```bibtex
> @article{PSCToolkit,
> title = {Parallel Sparse Computation Toolkit},
> journal = {Software Impacts},
> volume = {15},
> pages = {100463},
> year = {2023},
> issn = {2665-9638},
> doi = {https://doi.org/10.1016/j.simpa.2022.100463},
> url = {https://www.sciencedirect.com/science/article/pii/S2665963822001476},
> author = {Pasqua D’Ambra and Fabio Durastante and Salvatore Filippone},
> keywords = {Linear solvers, Algebraic preconditioners, HPC, GPU, Heterogeneous computing},
> }
> @techreport{RKToolbox,
>  title={A rational Krylov toolbox for MATLAB},
>  author={Berljafa, Mario and Elsworth, Steven and G{\"u}ttel, Stefan},
>  year={2014},
>  institution={Manchester Institute for Mathematical Sciences, University of Manchester},
>  url={http://guettel.com/rktoolbox/index.html},
> }
> @article{10.1145/2168773.2168781,
> author = {Niesen, Jitse and Wright, Will M.},
> title = {Algorithm 919: A Krylov Subspace Algorithm for Evaluating the ϕ-Functions Appearing in Exponential Integrators},
> year = {2012},
> issue_date = {April 2012},
> publisher = {Association for Computing Machinery},
> address = {New York, NY, USA},
> volume = {38},
> number = {3},
> issn = {0098-3500},
> url = {https://doi.org/10.1145/2168773.2168781},
> doi = {10.1145/2168773.2168781},
> journal = {ACM Trans. Math. Softw.},
> month = apr,
> articleno = {22},
> numpages = {19},
> keywords = {Krylov subspace methods, exponential integrators, matrix exponential}
> }
> ```

## Running the experiments

The code to run the experiments is contained in the :open_file_folder: matlabsrc folder. The folder also contains the various routines that implement the walk-based Laplacians. See the README inside the folder for more information.
