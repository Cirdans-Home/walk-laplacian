# Install

To install the PSBLAS3 and AMG4PSBLAS library here one could follow these steps from
the root directory of the repository.
```bash
cd psblas3
sh ../configure/psb_gnu.sh
make -j4 install
cd ..
cd amg4psblas
sh ../configure/amg_gnu.sh
make -j4 install
```
The configure options in the files `psb_gnu.sh` and `amg_gnu.sh` are set up for a
standard configuration with Ubuntu linux (v. 22.04). They can be modified to work
with different installations of the auxiliary and external libraries. Please
refer to the 
```bash
cd psblas3
./configure --help
```
and
```bash
cd amg4psblas
./configure --help
```
for further information on how to configure/compile/install the [PSCToolkit](https://psctoolkit.github.io/) on
different machines, some hints are also given in the configure folder.

