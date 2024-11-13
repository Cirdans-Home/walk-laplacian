PREFIX=$(cd ../install && pwd)

./configure --prefix=${PREFIX} \
	--with-psblas=${PREFIX}
