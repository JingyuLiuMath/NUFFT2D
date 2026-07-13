# Extern packages

## FFTW

Download `fftw-3.3.10.tar.gz` from [FFTW](https://fftw.org/), or use the following command:

``` bash
curl -O https://fftw.org/pub/fftw/fftw-3.3.10.tar.gz
```

Then use `tar -xvzf fftw-3.3.10.tar.gz` to extract it.

We need the single- and double-precision libraries of FFTW.

``` bash
cd fftw-3.3.10

./configure --prefix=<path/to/FFTWdouble> --enable-shared --enable-openmp
make -j
make install

./configure --enable-single --prefix=<path/to/FFTWsingle> --enable-shared --enable-openmp
make -j
make install
```

## FINUFFT

Download [FINUFFT](https://finufft.readthedocs.io/en/latest/index.html) by the following command:

``` bash
git clone https://github.com/flatironinstitute/finufft.git
```

Run the following commands:

``` bash
module load MATLAB
export PATH=$PATH:<path/to/FFTWdouble>/bin:<path/to/FFTWsingle>/bin
export LIBRARY_PATH=$LIBRARY_PATH:<path/to/FFTWdouble>/lib:<path/to/FFTWsingle>/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<path/to/FFTWdouble>/lib:<path/to/FFTWsingle>/lib
export CPATH=$CPATH:<path/to/FFTWdouble>/include:<path/to/FFTWsingle>/include

cd finufft
git checkout c6860a8d

cp make-platforms/make.inc.manylinux make.inc
make test -j
make matlab -j
```

Use

``` matlab
M = 1e5; x = 2*pi*rand(M,1); c = randn(M,1)+1i*randn(M,1); N = 2e5; f = finufft1d1(x,c,+1,1e-12,N);
```

to verify its success.

If you cannot get FINUFFT to compile, as a last resort you might find a precompiled binary for your platform under Assets for various [releases](https://github.com/flatironinstitute/finufft/releases).
