PREFIX=$(cd ../install && pwd)

./configure --prefix=${PREFIX} \
	--with-ccopt="-O3 -g -ggdb -fPIC -march=native -mtune=native" \
	--with-fcopt="-O3 -g -ggdb -frecursive -fPIC -march=native -mtune=native" \
	--with-cxxopt="-O3 -g -ggdb -fPIC -march=native -mtune=native" \
        --enable-cuda \
	--with-cudadir=$CUDA_HOME \
	--with-cudacc=86 \
	--enable-openmp

