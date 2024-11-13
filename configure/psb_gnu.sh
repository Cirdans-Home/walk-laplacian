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
