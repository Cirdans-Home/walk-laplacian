# Configure and Install PSCToolkit

The file `psb_gnu.sh` contains an example of how to configure the PSBLAS library on a machine having installed Ubuntu and the auxiliary libraries through system packages:
```bash
PREFIX=$(cd ../install && pwd)

./configure \
	--with-metisincdir=/opt/metis/include/ \
	--with-metislibdir=/opt/metis/lib/ \
	--with-ipk=4 --with-lpk=8 \
	--with-ccopt="-O3 -march=native"\
	--with-fcopt="-O3 -march=native"\
	--with-cxxopt="-O3 -march=native" \
	--with-amddir=/usr/include/suitesparse \
	--prefix=${PREFIX} \
	--enable-cuda \
	--with-cuda=${CUDA_HOME} \
	--with-cudacc=89
```
> [!WARNING]
> The examples in the paper make use of the library's GPU support, so they require an NVIDIA GPU with CUDA compute capabilities, the related libraries installed, and a compiler compatible with the CUDA version used. This usually add to the environment the `$CUDA_HOME` variable. If this is not the case, you have to manually tell to the configure where the CUDA libraries are installed.

After the configuration and installation of PSBLAS has been completed, you can move on to configuring AMG4PSBLAS, the library containing the preconditioner suite, using the `amg_gnu.sh` file:
```bash
PREFIX=$(cd ../install && pwd)

./configure --prefix=${PREFIX} \
        --with-psblas=${PREFIX} 
```
The test cases do not use external AMG4PSBLAS packages, so the required configuration is minimal.

## Installing on a cluster

The preferred environment for the PSCToolkit suite is a cluster where packages and libraries are typically made available through the use of environment modules. In particular, the experiments contained in the paper are performed using the Toeplitz cluster at the Green Data Center of the University of Pisa.

The file `environment_toeplitz.sh` contains the modules that are loaded to install the library.
```bash
module purge
module load gpu-gcc/12.2.0 gpu-metis/5.1.0-gcc-12.2.0 gpu-cuda/12.3.1-gcc-12.2.0 gpu-openmpi/4.1.6-cuda-12.3.1-gcc-12.2.0 gpu-openb
las/0.3.26-gcc-12.2.0
module list
```
Usually modules automatically set several environment variables, this simplifies the configuration phase, because the system is able to automatically find the auxiliary libraries.

See for example the `psb_gnu_toeplitz.sh` files:
```bash
PREFIX=$(cd ../install && pwd)

./configure --prefix=${PREFIX} \
        --with-ccopt="-O3 -g -ggdb -fPIC -march=native -mtune=native" \
        --with-fcopt="-O3 -g -ggdb -frecursive -fPIC -march=native -mtune=native" \
        --with-cxxopt="-O3 -g -ggdb -fPIC -march=native -mtune=native" \
        --enable-cuda \
        --with-cudadir=$CUDA_HOME \
        --with-cudacc=86 \
        --enable-openmp
```
and `amg_gnu_toeplitz.sh`:
```bash
PREFIX=$(cd ../install && pwd)

./configure --prefix=${PREFIX} \
        --with-psblas=${PREFIX} 
```

## Other machines
For installation on other machines or configurations, it is advisable to [:link: consult the documentation](https://psctoolkit.github.io/) and be guided by the help
```bash
./configure --help
```
in the two folders psblas3 and amg4psblas.
