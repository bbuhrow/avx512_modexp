/*
Copyright (c) 2021, Ben Buhrow
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.
*/

#ifndef _MONTY_H
#define _MONTY_H

#include "bigarith.h"

/* montgomery arithmetic operations */
typedef struct
{
	bignum *r;
	bignum *n;
    bignum *np;
	bignum *nhat;
    bignum *vnhat;
    bignum *rhat;
    bignum *rmask;
	bignum *one;
	bignum *mtmp1;
	bignum *mtmp2;
	bignum *mtmp3;
    bignum *mtmp4;
    base_t rho;    
} monty;

// montgomery arithmetic setup and conversion
void to_monty(monty *mdata, bignum * x);
monty * monty_alloc();
void monty_init(monty * in, bignum * n, int verbose);
void monty_free(monty *mdata);

// pointers to the current mul/sqr scanning technique
void(*mul_ptr)(monty *, bignum *, bignum *, bignum *, bignum *);
void(*sqr_ptr)(monty *, bignum *, bignum *, bignum *);

// montgomery multipliers for various scanning techniques
void mulmod_sos(monty *mdata, bignum * u, bignum * v, bignum * w, bignum * s);
void mulmod_cios(monty *mdata, bignum * u, bignum * v, bignum * w, bignum * s);
void mulmod_fios(monty *mdata, bignum * u, bignum * v, bignum * w, bignum * s);
void mulmod_fips(monty *mdata, bignum * u, bignum * v, bignum * w, bignum * s);
void mulmod_bps(monty* mdata, bignum* u, bignum* v, bignum* w, bignum* s);

// montgomery squaring that just calls mulmod for various scanning techniques
void sqrmod_sos_mul(monty *mdata, bignum * u, bignum * w, bignum * s);
void sqrmod_fips_mul(monty *mdata, bignum * u, bignum * w, bignum * s);
void sqrmod_fios_mul(monty *mdata, bignum * u, bignum * w, bignum * s);
void sqrmod_cios_mul(monty *mdata, bignum * u, bignum * w, bignum * s);
void sqrmod_bps_mul(monty* mdata, bignum* u, bignum* w, bignum* s);

// specialized montgomery squaring for various scanning techniques
void sqrmod_sos(monty *mdata, bignum * u, bignum * w, bignum * s);
void sqrmod_cios(monty* mdata, bignum* u, bignum* w, bignum* s);
void sqrmod_fios(monty *mdata, bignum * u, bignum * w, bignum * s);
void sqrmod_fips(monty* mdata, bignum* u, bignum* w, bignum* s);
void sqrmod_bps(monty* mdata, bignum* u, bignum* w, bignum* s);

#endif // _MONTY_H
