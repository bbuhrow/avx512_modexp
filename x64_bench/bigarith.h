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

#ifndef _BIG_ARITH_H
#define _BIG_ARITH_H

#include <stdint.h>
#include <stdlib.h>
#include "util.h"
#include "x64_arith.h"

#define DIGITBITS 64

#ifdef _MSC_VER
#define base_t unsigned long long
#define base_signed_t long long
#define MAXDIGIT 0xffffffffffffffff
#define HIBITMASK 0x8000000000000000
#define MAX_DEC_WORD 0x8AC7230489E80000ULL
#define DEC_DIGIT_PER_WORD 19
#define HEX_DIGIT_PER_WORD 16
#define HALFMASK 0xffffffff
#define HALFBITS 32
#else
#if DIGITBITS == 64

#define base_t uint64_t
#define base_signed_t int64_t
#define HALFBITS 32
#define HALFMASK 0xffffffff
#define MAXDIGIT 0xffffffffffffffff
#define HIBITMASK 0x8000000000000000
#define VECLEN 8
#define MAX_DEC_WORD 0x8AC7230489E80000ULL
#define DEC_DIGIT_PER_WORD 19
#define HEX_DIGIT_PER_WORD 16

#else

#define base_t uint32_t
#define base_signed_t int32_t
#define HALFBITS 16
#define HALFMASK 0xffff
#define MAXDIGIT 0xffffffff
#define HIBITMASK 0x80000000
#define VECLEN 16
#define MAX_DEC_WORD 0x3b9aca00
#define DEC_DIGIT_PER_WORD 9
#define HEX_DIGIT_PER_WORD 8

#endif
#endif

// supported: 
// any N divisible by 32 for MAXBITS < 512
// any N divisible by 128 MAXBITS >= 128
#ifndef MAXBITS
#define MAXBITS 512
#endif

#define NWORDS (MAXBITS / DIGITBITS)

typedef struct
{
	base_t *data;
	int size;
} bignum;


/* basic arithmetic operations: fixed allocation, variable sized, non-signed */
int zBits(bignum * n);
base_t spBits(base_t n);
void zSet1(bignum *dest, base_t value);
void zCopy(bignum * src, bignum * dest);
void zAdd(bignum * u, bignum * v, bignum * w);
void zShortAdd(bignum * u, base_t v, bignum * w);
int zSub(bignum * u, bignum * v, bignum * w);
void zShortSub(bignum * u, base_t v, bignum * w);
int zCompare(bignum * u, bignum * v);
int zCompare1(bignum * u, base_t v);
base_t zShortDiv(bignum * u, base_t v, bignum * q);
void zDiv(bignum * u, bignum * v, bignum * q, bignum * r);
int shortCompare(base_t p[2], base_t t[2]);
int shortSubtract(base_t u[2], base_t v[2], base_t w[2]);
void zMul(bignum * u, bignum * v, bignum * w);
void zMult(bignum * u, bignum * v, bignum * w, bignum *tmp);
void zModMul(bignum * u, bignum * v, bignum * n, bignum * w);
void zShortMul(bignum * u, base_t v, bignum * w);
void zSqr(bignum * x, bignum * w);
void zShiftLeft(bignum * a, bignum * b, int x);
void zShiftLeft_1(bignum * a, bignum * b);
void zShiftRight(bignum * a, bignum * b, int x);
void zShiftRight_1(bignum * a, bignum * b);
void spAdd(base_t u, base_t v, base_t *sum, base_t *carry);
void spAdd3(base_t u, base_t v, base_t w, base_t *sum, base_t *carry);
void spSub3(base_t u, base_t v, base_t w, base_t *sub, base_t *borrow);
void spSub(base_t u, base_t v, base_t *sub, base_t *borrow);
base_t spDivide(base_t *q, base_t *r, base_t u[2], base_t v);
void spMultiply(base_t u, base_t v, base_t *product, base_t *carry);
void spMulAdd(base_t u, base_t v, base_t w, base_t t, base_t *lower, base_t *carry);
void spMulMod(base_t u, base_t v, base_t m, base_t *w);
void sp2big(base_t src, bignum * dest);
void zClear(bignum * n);
void zClearFull(bignum * n);
void zClamp(bignum * n);
void zPrint(bignum *n);
int ndigits_1(base_t n);
bignum * zInit(void);
void zFree(bignum *n);
void xGCD(bignum *a, bignum *b, bignum *x, bignum *y, bignum *g);
int zBinGCD(bignum *u, bignum *v, bignum *w);
int zLEGCD(bignum *u, bignum *v, bignum *w);
base_t spGCD(base_t x, base_t y);
void zModMuls(bignum * u, bignum * v, bignum * n, bignum * w, bignum *s1, bignum *s2);
void zModExp(bignum *d, bignum *b, bignum *e, bignum *m);

void str2hexz(char in[], bignum * u);
void zDec2Hex(bignum * u, bignum * v);
char *z2decstr(bignum * n);
void zHex2Dec(bignum * u, bignum * v);

#endif // _BIGARITH_H

