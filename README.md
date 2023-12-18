# avx512_modexp
A test library for computing modular exponentiation in parallel using AVX-512 vector arithmetic

Verified to work with gcc 7.3.0, gcc 11.1.0 and icc 18.0.3.

build with (required)
make SKYLAKEX=1

optionally add this to build line to change the length of test numbers (N needs to be a multiple of 128)
MAXBITS=N

optionally add this to build line to change the compiler to gcc-7.3.0 from the default icc
COMPILER=gcc730

optionally add this to build line to change the compiler to gcc-11.1.0 from the default icc
COMPILER=gcc11

optionally add this to build line to use the double precision FMA arithmetic instead of 32-bit integer arithmetic.
If this is specified then MAXBITS must be multiples of 208
BASE51=1


Run the executable with 2 arguments: number of threads and whether or not to verify all of the results using GMP.
For example, to run with 4 threads and skip verification:
./avx512_modexp 4 0
