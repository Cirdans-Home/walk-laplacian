# Data

The script in this folder can be used to download the test matrices/complex
network from the [SuiteSparse collection](https://sparse.tamu.edu). The file in the `mm` subfolder 
performs also the conversion to Matrix Market format of the largest connected
component of the graph. These files are used for the Fortran testers.

## Cite data source

If you use this Collection as a data set, please cite the original ACM-TOMS paper:
```bibtex
@article{10.1145/2049662.2049663,
author = {Davis, Timothy A. and Hu, Yifan},
title = {The University of Florida Sparse Matrix Collection},
year = {2011},
issue_date = {November 2011},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
volume = {38},
number = {1},
issn = {0098-3500},
url = {https://doi.org/10.1145/2049662.2049663},
doi = {10.1145/2049662.2049663},
abstract = {We describe the University of Florida Sparse Matrix Collection, a large and actively growing set of sparse matrices that arise in real applications. The Collection is widely used by the numerical linear algebra community for the development and performance evaluation of sparse matrix algorithms. It allows for robust and repeatable experiments: robust because performance results with artificially generated matrices can be misleading, and repeatable because matrices are curated and made publicly available in many formats. Its matrices cover a wide spectrum of domains, include those arising from problems with underlying 2D or 3D geometry (as structural engineering, computational fluid dynamics, model reduction, electromagnetics, semiconductor devices, thermodynamics, materials, acoustics, computer graphics/vision, robotics/kinematics, and other discretizations) and those that typically do not have such geometry (optimization, circuit simulation, economic and financial modeling, theoretical and quantum chemistry, chemical process simulation, mathematics and statistics, power networks, and other networks and graphs). We provide software for accessing and managing the Collection, from MATLAB™, Mathematica™, Fortran, and C, as well as an online search capability. Graph visualization of the matrices is provided, and a new multilevel coarsening scheme is proposed to facilitate this task.},
journal = {ACM Trans. Math. Softw.},
month = {dec},
articleno = {1},
numpages = {25},
keywords = {Graph drawing, performance evaluation, multilevel algorithms, sparse matrices}
}
```
