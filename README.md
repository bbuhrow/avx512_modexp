# avx512_modexp
A test library for computing modular exponentiation in parallel using AVX-512 vector arithmetic

Verified to work with gcc 7.3.0 and icc 18.0.3.

build with (required)
make SKYLAKEX=1

optionally add this to build line to change the length of test numbers (N needs to be a multiple of 128)
MAXBITS=N

optionally add this to build line to change the compiler to gcc-7.3.0 from the default icc (I've tested only gcc730 and icc)
COMPILER=gcc730
