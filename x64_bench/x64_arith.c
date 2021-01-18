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
#include "x64_arith.h"

__inline void spAdd(uint64_t u, uint64_t v, uint64_t *sum, uint64_t *carry)
{
    uint64_t s, c;

    s = v;
    c = 0;

    __asm__("movq %2, %%rax		\n\t"
        "addq %%rax, %3		\n\t"
        "adcq $0, %4		\n\t"
        : "=r"(s), "=r"(c)
        : "r"(u), "0"(s), "1"(c)
        : "%rax", "memory", "cc");

    *sum = s;
    *carry = c;

    return;
}

__inline void spAdd3(uint64_t u, uint64_t v, uint64_t w, uint64_t *sum, uint64_t *carry)
{
    uint64_t s, c;

    s = v;
    c = 0;

    __asm__("movq %2, %%rax		\n\t"
        "addq %3, %%rax		\n\t"
        "adcq $0, %5		\n\t"
        "addq %%rax, %4		\n\t"
        "adcq $0, %5		\n\t"
        : "=r"(s), "=r"(c)
        : "r"(u), "r"(w), "0"(s), "1"(c)
        : "%rax", "memory", "cc");

    *sum = s;
    *carry = c;

    return;
}

__inline void spSub3(uint64_t u, uint64_t v, uint64_t w, uint64_t *sub, uint64_t *borrow)
{
    uint64_t s, b;

    s = v;
    b = 0;

    __asm__("movq %2, %%rax		\n\t"
        "subq %4, %%rax		\n\t"
        "adcq $0, %5		\n\t"
        "subq %3, %%rax		\n\t"
        "adcq $0, %5		\n\t"
        "movq %%rax, %4		\n\t"
        : "=r"(s), "=r"(b)
        : "r"(u), "r"(w), "0"(s), "1"(b)
        : "%rax", "memory", "cc");

    *sub = s;
    *borrow = b;

    return;
}

__inline void spSub(uint64_t u, uint64_t v, uint64_t *sub, uint64_t *borrow)
{
    uint64_t s, b;

    s = v;
    b = 0;

    __asm__("movq %2, %%rax		\n\t"
        "subq %3, %%rax		\n\t"
        "adcq $0, %4		\n\t"
        "movq %%rax, %3		\n\t"
        : "=r"(s), "=r"(b)
        : "r"(u), "0"(s), "1"(b)
        : "%rax", "memory", "cc");

    *sub = s;
    *borrow = b;

    return;
}

__inline uint64_t spDivide(uint64_t *q, uint64_t *r, uint64_t u[2], uint64_t v)
{
    *r = u[1];
    *q = u[0];
    __asm__("divq %4"
        : "=a"(*q), "=d"(*r)
        : "1"(*r), "0"(*q), "r"(v));

    return 0;
}

__inline void spMultiply(uint64_t u, uint64_t v, uint64_t *product, uint64_t *carry)
{
    *product = v;
    *carry = u;

    __asm__("movq %2, %%rax	\n\t"
        "mulq %3	\n\t"
        "movq %%rax, %0		\n\t"
        "movq %%rdx, %1		\n\t"
        : "=r"(*product), "=r"(*carry)
        : "1"(*carry), "0"(*product)
        : "%rax", "%rdx", "memory", "cc");

    return;
}

__inline uint64_t spDiv(uint64_t *q, uint64_t *r, uint64_t u1, uint64_t u0, uint64_t v)
{
    *r = u1;
    *q = u0;
    __asm__("divq %4"
        : "=a"(*q), "=d"(*r)
        : "1"(*r), "0"(*q), "r"(v));

    return 0;
}

__inline uint64_t spMod(uint64_t u1, uint64_t u0, uint64_t v)
{
    __asm__("divq %4"
        : "=a"(u0), "=d"(u1)
        : "1"(u1), "0"(u0), "r"(v));

    return u1;
}

__inline void spMul(uint64_t u, uint64_t v, uint64_t *product, uint64_t *carry)
{
    *product = v;
    *carry = u;

    __asm__("movq %2, %%rax	\n\t"
        "mulq %3	\n\t"
        "movq %%rax, %0		\n\t"
        "movq %%rdx, %1		\n\t"
        : "=r"(*product), "=r"(*carry)
        : "0"(*product), "1"(*carry)
        : "rax", "rdx", "cc");

    return;
}

__inline void spMulAdd1(uint64_t u, uint64_t v, uint64_t w,
    uint64_t *product, uint64_t *carry)
{
    *product = v;
    *carry = u;

    __asm__("movq %2, %%rax	\n\t"
        "mulq %3	\n\t"
        "addq %4, %%rax \n\t"
        "adcq $0, %%rdx \n\t"
        "movq %%rax, %0		\n\t"
        "movq %%rdx, %1		\n\t"
        : "=r"(*product), "=r"(*carry)
        : "1"(*carry), "0"(*product), "r"(w)
        : "rax", "rdx", "cc");

    return;
}

__inline void spMulAdd2(uint64_t u, uint64_t v, uint64_t w,
    uint64_t c, uint64_t *product, uint64_t *carry)
{
    *product = v;
    *carry = u;

    __asm__("movq %2, %%rax	\n\t"
        "mulq %3	\n\t"
        "addq %4, %%rax \n\t"
        "adcq $0, %%rdx \n\t"
        "addq %5, %%rax \n\t"
        "adcq $0, %%rdx \n\t"
        "movq %%rax, %0		\n\t"
        "movq %%rdx, %1		\n\t"
        : "=r"(*product), "=r"(*carry)
        : "1"(*carry), "0"(*product), "r"(w), "r"(c)
        : "rax", "rdx", "cc");

    return;
}

__inline void spMulAdd2x(uint64_t u, uint64_t v, uint64_t w,
    uint64_t c, uint64_t *product, uint64_t *carry)
{
    *product = v;
    *carry = u;

    // maximum in all inputs won't overflow outputs:
    // 0xffffffffffffffff ^ 2 + 2 * 0xffffffffffffffff = 0xffffffffffffffffffffffffffffffff

    __asm__("movq %2, %%rdx	\n\t"
        "addq %5, %4 \n\t"          /* add current output to previous carry */
        "mulx %3, %0, %1	\n\t"   /* multiply */
        "adcq $0, %1 \n\t"          /* add carry into himul result */
        "addq %4, %0 \n\t"          /* lowmul + current output + previous carry, store into current output */
        "adcq $0, %1 \n\t"          /* carry prop into himul */
        : "=r"(*product), "=r"(*carry)
        : "1"(*carry), "0"(*product), "r"(w), "r"(c)
        : "r10", "rdx", "r11", "r12", "cc");

    return;
}

__inline void mpSub(uint64_t * u, uint64_t * n, uint64_t * w, int sz)
{
    int i;
    uint64_t b, d;

    b = 0;
    for (i = 0; i < sz; i++)
    {
        spSub3(u[i], n[i], b, &w[i], &b);
    }

    if (b)
        spSub(u[i], b, &w[i], &b);

    return;
}

__inline void mpSub1(uint64_t * u, uint64_t n, uint64_t * w, int sz)
{
    int i = 0;
    uint64_t b;

    b = 0;
    spSub3(u[i], n, b, &w[i], &b);
    i++;
    while (i < sz)
    {
        spSub(u[i], b, &w[i], &b);
        i++;
    }

    if (b)
        spSub(u[i], b, &w[i], &b);

    return;
}

__inline void spMulAddc(uint64_t u, uint64_t v, uint64_t * w)
{
    // for use with product scanning approach...
    // multiply u*v.
    // add result into w[0] and w[1] and carry propagate once.

    __asm__("movq %0, %%rax	\n\t"
        "mulq %1	\n\t"
        "movq 16(%2), %%r10	\n\t"
        "addq 0(%2), %%rax \n\t"
        "adcq 8(%2), %%rdx \n\t"
        "adcq $0, %%r10 \n\t"
        "movq %%rax, 0(%2)		\n\t"
        "movq %%rdx, 8(%2)		\n\t"
        "movq %%r10, 16(%2)		\n\t"
        :
    : "r"(u), "r"(v), "r"(w)
        : "rax", "rdx", "r10", "cc", "memory");

    return;
}

__inline void spMul2Acc(uint64_t u, uint64_t v, uint64_t n, uint64_t s, uint64_t * w)
{
    // for use with product scanning approach...
    // multiply u*v.
    // add result into w[0] and w[1] and carry propagate once.
    // multiply n*s.
    // add result into w[0] and w[1] and carry propagate once.

    __asm__("movq %0, %%rax	\n\t"
        "mulq %1	\n\t"
        "movq 16(%4), %%r10	\n\t"
        "addq 0(%4), %%rax \n\t"
        "movq %%rax, %%r11 \n\t"
        "adcq 8(%4), %%rdx \n\t"
        "movq %%rdx, %%r12 \n\t"
        "adcq $0, %%r10 \n\t"
        "movq %2, %%rax	\n\t"
        "mulq %3	\n\t"
        "addq %%r11, %%rax \n\t"
        "adcq %%r12, %%rdx \n\t"
        "adcq $0, %%r10 \n\t"
        "movq %%rax, 0(%4)		\n\t"
        "movq %%rdx, 8(%4)		\n\t"
        "movq %%r10, 16(%4)		\n\t"
        :
    : "r"(u), "r"(v), "r"(n), "r"(s), "r"(w)
        : "rax", "rdx", "r10", "r11", "r12", "cc", "memory");

    return;
}

__inline void spMulAddcr(uint64_t u, uint64_t v, uint64_t * w)
{
    // for use with product scanning approach...
    // multiply u*v.
    // add result into w[0] and w[1] and carry propagate once.
    // final output rotation.

    __asm__("movq %0, %%rax	\n\t"
        "mulq %1	\n\t"
        "movq 16(%2), %%r10	\n\t"
        "addq 0(%2), %%rax \n\t"
        "adcq 8(%2), %%rdx \n\t"
        "adcq $0, %%r10 \n\t"
        "xorq %%rax, %%rax \n\t"
        "movq %%rdx, 0(%2)		\n\t"
        "movq %%r10, 8(%2)		\n\t"
        "movq %%rax, 16(%2)		\n\t"
        :
    : "r"(u), "r"(v), "r"(w)
        : "rax", "rdx", "r10", "cc", "memory");

    return;
}

__inline void spMulDblAdd_1(uint64_t u, uint64_t v, uint64_t carryin, uint64_t * w, uint64_t *carryout)
{
    // for use with sos squaring approach...
    // multiply u*v and add carryin to the 2nd result word.
    // add result twice into w[0] and w[1] and return any further carryout.

    __asm__("movq %3, %%rax	\n\t"
        "mulq %4	\n\t"
        "xorq %%r10, %%r10	\n\t"
        "addq %%rax, %%rax \n\t"
        "adcq %%rdx, %%rdx \n\t"
        "adcq $0, %%r10 \n\t"
        "addq %5, %%rax \n\t"
        "adcq %6, %%rdx \n\t"
        "adcq $0, %%r10 \n\t"
        "adcq %7, %%rdx \n\t"
        "adcq $0, %%r10 \n\t"
        "movq %%rax, %0		\n\t"
        "movq %%rdx, %1		\n\t"
        "movq %%r10, %2		\n\t"
        : "=r"(w[0]), "=r"(w[1]), "=r"(*carryout)
        : "r"(u), "r"(v), "0"(w[0]), "1"(w[1]), "r"(carryin)
        : "rax", "rdx", "r10", "cc");

    return;
}

__inline void spMulDblAdd_2(uint64_t u, uint64_t v, uint64_t carryin, uint64_t * w, uint64_t *carryout)
{
    // for use with sos squaring approach...
    // multiply u*v and add carryin to the 2nd result word.
    // add result twice into w[0] and w[1] and return any further carryout.
    // same approach as _1, except instead of add and adc we use shldq/shl

    __asm__("movq %3, %%rax	\n\t"
        "mulq %4	\n\t"
        "xorq %%r10, %%r10	\n\t"
        "shldq $1, %%rax, %%rdx \n\t"
        "adcq $0, %%r10 \n\t"
        "shlq $1, %%rax \n\t"
        "addq %5, %%rax \n\t"
        "adcq %6, %%rdx \n\t"
        "adcq $0, %%r10 \n\t"
        "adcq %7, %%rdx \n\t"
        "adcq $0, %%r10 \n\t"
        "movq %%rax, %0		\n\t"
        "movq %%rdx, %1		\n\t"
        "movq %%r10, %2		\n\t"
        : "=r"(w[0]), "=r"(w[1]), "=r"(*carryout)
        : "r"(u), "r"(v), "0"(w[0]), "1"(w[1]), "r"(carryin)
        : "rax", "rdx", "r10", "cc");

    return;
}

__inline void spMulDblAdd_3(uint64_t u, uint64_t v, uint64_t * w)
{
    // for use with fips squaring approach...
    // multiply u*v.
    // add result twice into w[0], w[1], and w[2].

    __asm__("movq %3, %%rax	\n\t"
        "mulq %4	\n\t"
        "movq %7, %%r10	\n\t"
        "addq %%rax, %%rax \n\t"
        "adcq %%rdx, %%rdx \n\t"
        "adcq $0, %%r10 \n\t"
        "addq %5, %%rax \n\t"
        "adcq %6, %%rdx \n\t"
        "adcq $0, %%r10 \n\t"
        "movq %%rax, %0		\n\t"
        "movq %%rdx, %1		\n\t"
        "movq %%r10, %2		\n\t"
        : "=r"(w[0]), "=r"(w[1]), "=r"(w[2])
        : "r"(u), "r"(v), "0"(w[0]), "1"(w[1]), "r"(w[2])
        : "rax", "rdx", "r10", "cc");

    return;
}

__inline void spSqrMulAcc(uint64_t u, uint64_t v, uint64_t n, uint64_t s, uint64_t * w)
{
    // for use with fips squaring approach on cross-terms...
    // multiply u*v.
    // add result twice into w[0], w[1], and w[2].
    // multiply n*s.
    // add result once into w[0], w[1], and w[2].

    __asm__("movq %3, %%rax	\n\t"
        "mulq %4	\n\t"
        "movq %7, %%r10	\n\t"
        "addq %%rax, %%rax \n\t"
        "adcq %%rdx, %%rdx \n\t"
        "adcq $0, %%r10 \n\t"
        "addq %5, %%rax \n\t"
        "movq %%rax, %%r11 \n\t"
        "adcq %6, %%rdx \n\t"
        "movq %%rdx, %%r12 \n\t"
        "adcq $0, %%r10 \n\t"
        "movq %8, %%rax	\n\t"
        "mulq %9	\n\t"
        "addq %%r11, %%rax \n\t"
        "adcq %%r12, %%rdx \n\t"
        "adcq $0, %%r10 \n\t"
        "movq %%rax, %0		\n\t"
        "movq %%rdx, %1		\n\t"
        "movq %%r10, %2		\n\t"
        : "=r"(w[0]), "=r"(w[1]), "=r"(w[2])
        : "r"(u), "r"(v), "0"(w[0]), "1"(w[1]), "r"(w[2]), "r"(n), "r"(s)
        : "rax", "rdx", "r10", "r11", "r12", "cc");

    return;
}

void mpAdd1b(uint64_t * u, uint64_t n, uint64_t * w, int sz)
{
    // assume u and w point to the same thing, so we
    // can stop as soon as there is no carry.
    int i = 0;
    uint64_t c = 0;

    spAdd3(u[i], n, c, &w[i], &c);
    i++;
    while ((i < sz) && (c > 0))
    {
        spAdd(u[i], c, &w[i], &c);
        i++;
    }

    if (c)
        spAdd(u[i], c, &w[i], &c);

    return;
}

void mpAdd1(uint64_t * u, uint64_t n, uint64_t * w, int sz)
{
    int i = 0;
    uint64_t c;

    c = 0;
    spAdd3(u[i], n, c, &w[i], &c);
    i++;
    while (i < sz)
    {
        spAdd(u[i], c, &w[i], &c);
        i++;
    }

    if (c)
        spAdd(u[i], c, &w[i], &c);

    return;
}

void mpAdd(uint64_t * u, uint64_t * v, uint64_t * w, int sz)
{
    int i = 0;
    uint64_t c;

    c = 0;
    for (i = 0; i < sz; i++)
    {
        spAdd3(u[i], v[i], c, &w[i], &c);
    }
    w[i] = c;

    return;
}

