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

#ifndef _X64_ARITH_H
#define _X64_ARITH_H

// this file declares special routines for low-level x64 arithmetic
// used as subroutines for modular multiplication.
#include <stdint.h>




void spAdd(uint64_t u, uint64_t v, uint64_t *sum, uint64_t *carry);
void spAdd3(uint64_t u, uint64_t v, uint64_t w, uint64_t *sum, uint64_t *carry);
void spSub3(uint64_t u, uint64_t v, uint64_t w, uint64_t *sub, uint64_t *borrow);
void spSub(uint64_t u, uint64_t v, uint64_t *sub, uint64_t *borrow);
uint64_t spDivide(uint64_t *q, uint64_t *r, uint64_t u[2], uint64_t v);
void spMultiply(uint64_t u, uint64_t v, uint64_t *product, uint64_t *carry);
uint64_t spDiv(uint64_t *q, uint64_t *r, uint64_t u1, uint64_t u0, uint64_t v);
uint64_t spMod(uint64_t u1, uint64_t u0, uint64_t v);
void spMul(uint64_t u, uint64_t v, uint64_t *product, uint64_t *carry);
void spMulAdd1(uint64_t u, uint64_t v, uint64_t w,
    uint64_t *product, uint64_t *carry);
void spMulAdd2(uint64_t u, uint64_t v, uint64_t w,
    uint64_t c, uint64_t *product, uint64_t *carry);
void spMulAdd2x(uint64_t u, uint64_t v, uint64_t w,
    uint64_t c, uint64_t *product, uint64_t *carry);
void spMulAddc(uint64_t u, uint64_t v, uint64_t * w);
void spSqrMulAcc(uint64_t u, uint64_t v, uint64_t n, uint64_t s, uint64_t * w);
void mpSub(uint64_t * u, uint64_t * n, uint64_t * w, int sz);
void mpSub1(uint64_t * u, uint64_t n, uint64_t * w, int sz);
void mpAdd1(uint64_t * u, uint64_t n, uint64_t * w, int sz);
void spMul2Acc(uint64_t u, uint64_t v, uint64_t n, uint64_t s, uint64_t * w);
void spMulAddcr(uint64_t u, uint64_t v, uint64_t * w);
void spMulDblAdd_1(uint64_t u, uint64_t v, uint64_t carryin, uint64_t * w, uint64_t *carryout);
void spMulDblAdd_2(uint64_t u, uint64_t v, uint64_t carryin, uint64_t * w, uint64_t *carryout);
void spMulDblAdd_3(uint64_t u, uint64_t v, uint64_t * w);
void mpAdd1b(uint64_t * u, uint64_t n, uint64_t * w, int sz);


#endif
