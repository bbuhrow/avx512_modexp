/*
Copyright (c) 2014, Ben Buhrow
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


Copyright (c) 2018 by The Mayo Clinic, though its Special Purpose
 Processor Development Group (SPPDG). All Rights Reserved Worldwide.
 Licensed under the Apache License, Version 2.0 (the "License"); you may
 not use this file except in compliance with the License. You may obtain
 a copy of the License at http://www.apache.org/licenses/LICENSE-2.0.
Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied,
           including conditions of title, non-infringement, merchantability,
          or fitness for a particular purpose
 See the License for the specific language governing permissions and
 limitations under the License.
This file is a snapshot of a work in progress, originated by Mayo
 Clinic SPPDG.
*/

#include "vecarith.h"

#define MUL_ACCUM_2X(a0, b0, s0, n0, hiword) \
    prod1 = _mm512_mul_epu32(a0, b0);                          \
    prod2 = _mm512_mul_epu32(s0, n0);                          \
                                                                 \
    te0 = _mm512_add_epi64(te0, prod1);                        \
    scarry_1 = _mm512_cmplt_epu64_mask(te0, prod1);           \
    te0 = _mm512_add_epi64(te0, prod2);                        \
    scarry_2 = _mm512_cmplt_epu64_mask(te0, prod2);           \
    te1 = _mm512_mask_add_epi64(te1, scarry_1, hiword, te1);    \
    te1 = _mm512_mask_add_epi64(te1, scarry_2, hiword, te1);    \
                                                                 \
    a0 = _mm512_shuffle_epi32(a0, 0xB1);                            \
    b0 = _mm512_shuffle_epi32(b0, 0xB1);                            \
    s0 = _mm512_shuffle_epi32(s0, 0xB1);                            \
    n0 = _mm512_shuffle_epi32(n0, 0xB1);                            \
                                                                 \
    prod1 = _mm512_mul_epu32(a0, b0);                          \
    prod2 = _mm512_mul_epu32(s0, n0);                          \
                                                                 \
    to0 = _mm512_add_epi64(to0, prod1);                        \
    scarry_1 = _mm512_cmplt_epu64_mask(to0, prod1);           \
    to0 = _mm512_add_epi64(to0, prod2);                        \
    scarry_2 = _mm512_cmplt_epu64_mask(to0, prod2);           \
    to1 = _mm512_mask_add_epi64(to1, scarry_1, hiword, to1);    \
    to1 = _mm512_mask_add_epi64(to1, scarry_2, hiword, to1);

#define MUL_ACCUM_4X(a0, b0, b1, b2, b3, hiword) \
    prod1_e = _mm512_mul_epu32(a0, b0);                          \
    prod2_e = _mm512_mul_epu32(a0, b1);                          \
    prod3_e = _mm512_mul_epu32(a0, b2);                          \
    prod4_e = _mm512_mul_epu32(a0, b3);                          \
                                                                 \
    te0 = _mm512_add_epi64(te0, prod1_e);                        \
    te2 = _mm512_add_epi64(te2, prod2_e);                        \
    te4 = _mm512_add_epi64(te4, prod3_e);                        \
    te6 = _mm512_add_epi64(te6, prod4_e);                        \
    scarry_e1 = _mm512_cmplt_epu64_mask(te0, prod1_e);           \
    scarry_e2 = _mm512_cmplt_epu64_mask(te2, prod2_e);           \
    scarry_e3 = _mm512_cmplt_epu64_mask(te4, prod3_e);           \
    scarry_e4 = _mm512_cmplt_epu64_mask(te6, prod4_e);           \
    te1 = _mm512_mask_add_epi64(te1, scarry_e1, hiword, te1);    \
    te3 = _mm512_mask_add_epi64(te3, scarry_e2, hiword, te3);    \
    te5 = _mm512_mask_add_epi64(te5, scarry_e3, hiword, te5);    \
    te7 = _mm512_mask_add_epi64(te7, scarry_e4, hiword, te7);    \
                                                                 \
    prod1_e = _mm512_shuffle_epi32(b0, 0xB1);                    \
    prod2_e = _mm512_shuffle_epi32(b1, 0xB1);                    \
    prod3_e = _mm512_shuffle_epi32(b2, 0xB1);                    \
    prod4_e = _mm512_shuffle_epi32(b3, 0xB1);                    \
    prod5_e = _mm512_shuffle_epi32(a0, 0xB1);                    \
                                                                 \
    prod1_e = _mm512_mul_epu32(prod5_e, prod1_e);                \
    prod2_e = _mm512_mul_epu32(prod5_e, prod2_e);                \
    prod3_e = _mm512_mul_epu32(prod5_e, prod3_e);                \
    prod4_e = _mm512_mul_epu32(prod5_e, prod4_e);                \
                                                                 \
    to0 = _mm512_add_epi64(to0, prod1_e);                        \
    to2 = _mm512_add_epi64(to2, prod2_e);                        \
    to4 = _mm512_add_epi64(to4, prod3_e);                        \
    to6 = _mm512_add_epi64(to6, prod4_e);                        \
    scarry_e1 = _mm512_cmplt_epu64_mask(to0, prod1_e);           \
    scarry_e2 = _mm512_cmplt_epu64_mask(to2, prod2_e);           \
    scarry_e3 = _mm512_cmplt_epu64_mask(to4, prod3_e);           \
    scarry_e4 = _mm512_cmplt_epu64_mask(to6, prod4_e);           \
    to1 = _mm512_mask_add_epi64(to1, scarry_e1, hiword, to1);    \
    to3 = _mm512_mask_add_epi64(to3, scarry_e2, hiword, to3);    \
    to5 = _mm512_mask_add_epi64(to5, scarry_e3, hiword, to5);    \
    to7 = _mm512_mask_add_epi64(to7, scarry_e4, hiword, to7);

#define ACCUM_EO_PROD(sum_e, sum_o, carry_e, carry_o) \
	sum_e = _mm512_add_epi64(sum_e, prod1_e);	\
	sum_o = _mm512_add_epi64(sum_o, prod1_o);	\
	scarry_e1 = _mm512_cmplt_epu64_mask(sum_e, prod1_e);	\
	scarry_o1 = _mm512_cmplt_epu64_mask(sum_o, prod1_o);	\
	carry_e = _mm512_mask_add_epi64(carry_e, scarry_e1, hiword, carry_e);	\
	carry_o = _mm512_mask_add_epi64(carry_o, scarry_o1, hiword, carry_o);

#define ACCUM_EO_PROD2(sum_e, sum_o, carry_e, carry_o, in_e, in_o) \
	sum_e = _mm512_add_epi64(sum_e, in_e);	\
    sum_o = _mm512_add_epi64(sum_o, in_o);	\
	scarry_e1 = _mm512_cmplt_epu64_mask(sum_e, in_e);	\
    scarry_o1 = _mm512_cmplt_epu64_mask(sum_o, in_o);	\
	carry_e = _mm512_mask_add_epi64(carry_e, scarry_e1, hiword, carry_e); \
    carry_o = _mm512_mask_add_epi64(carry_o, scarry_o1, hiword, carry_o);

#define ACCUM_DOUBLED_EO_PROD(sum_e, sum_o, carry_e, carry_o) \
	sum_e = _mm512_add_epi64(sum_e, prod1_e);	\
	sum_o = _mm512_add_epi64(sum_o, prod1_o);	\
	scarry_e1 = _mm512_cmplt_epu64_mask(sum_e, prod1_e);	\
	scarry_o1 = _mm512_cmplt_epu64_mask(sum_o, prod1_o);	\
	carry_e = _mm512_mask_add_epi64(carry_e, scarry_e1, hiword2, carry_e);	\
	carry_o = _mm512_mask_add_epi64(carry_o, scarry_o1, hiword2, carry_o);

#define BLOCKWORDS 4
#define NBLOCKS (NWORDS / BLOCKWORDS)


void vecmulmod_bps(bignum *a, bignum *b, bignum *c, bignum *n, bignum *s, monty *mdata)
{
    int i, j, k;

    __m512i a0, a1, a2, a3;
    __m512i b0, b1, b2, b3;
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;
    __m512i to0, to1, to2, to3, to4, to5, to6, to7;

    __m512i acc_e0;
    __m512i acc_o0;

    __m512i acc_e1;
    __m512i acc_o1;

    __m512i nhatvec_e = _mm512_load_epi32(mdata->vrho);
    __m512i nhatvec_o = _mm512_shuffle_epi32(nhatvec_e, 0xB1);;

    __m512i prod1_e;
    __m512i prod1_o;

    __m512i hiword = _mm512_set1_epi64(0x000000100000000);
    __m512i zero = _mm512_set1_epi64(0);

    __mmask8 scarry_e1 = 0;
    __mmask8 scarry_o1 = 0;
    __mmask16 scarry2;
    __mmask16 scarry;

    // zero the accumulator
    acc_e0 = acc_o0 = acc_e1 = acc_o1 = zero;


    // first half mul
    for (i = 0; i < NBLOCKS; i++)
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;
        to0 = to1 = to2 = to3 = to4 = to5 = to6 = to7 = zero;

        for (j = i; j > 0; j--)
        {
            a0 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);

            __m512i prod1_e, prod2_e, prod3_e, prod4_e, prod5_e;
            __mmask8 scarry_e1, scarry_e2, scarry_e3, scarry_e4;

            // a0 terms
            MUL_ACCUM_4X(a0, b0, b1, b2, b3, hiword);

            // a1 terms
            b0 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            MUL_ACCUM_4X(a1, b1, b2, b3, b0, hiword);

            // a2 terms
            b1 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            MUL_ACCUM_4X(a2, b2, b3, b0, b1, hiword);

            // a3 terms
            b2 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);
            MUL_ACCUM_4X(a3, b3, b0, b1, b2, hiword);
        }

        // finish each triangular shaped column sum
        a0 = _mm512_load_epi32(a->data + (i * BLOCKWORDS + 0) * VECLEN);
        a1 = _mm512_load_epi32(a->data + (i * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi32(a->data + (i * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi32(a->data + (i * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi32(b->data + 0 * VECLEN);
        b1 = _mm512_load_epi32(b->data + 1 * VECLEN);
        b2 = _mm512_load_epi32(b->data + 2 * VECLEN);
        b3 = _mm512_load_epi32(b->data + 3 * VECLEN);

        // ======
        _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        // ======
        _mm512_mul_eo64_epi32(a1, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te2, to2, te3, to3);

        _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te2, to2, te3, to3);

        // ======
        _mm512_mul_eo64_epi32(a2, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te4, to4, te5, to5);

        _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te4, to4, te5, to5);

        _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te4, to4, te5, to5);

        // ======
        _mm512_mul_eo64_epi32(a3, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te6, to6, te7, to7);

        _mm512_mul_eo64_epi32(a2, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te6, to6, te7, to7);

        _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te6, to6, te7, to7);

        _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te6, to6, te7, to7);

        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        for (j = i; j > 0; j--)
        {
            a0 = _mm512_load_epi32(s->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi32(s->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi32(s->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi32(s->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);

            __m512i prod1_e, prod2_e, prod3_e, prod4_e, prod5_e;
            __mmask8 scarry_e1, scarry_e2, scarry_e3, scarry_e4;

            // a0 terms
            MUL_ACCUM_4X(a0, b0, b1, b2, b3, hiword);

            // a1 terms
            b0 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            MUL_ACCUM_4X(a1, b1, b2, b3, b0, hiword);

            // a2 terms
            b1 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            MUL_ACCUM_4X(a2, b2, b3, b0, b1, hiword);

            // a3 terms
            b2 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);
            MUL_ACCUM_4X(a3, b3, b0, b1, b2, hiword);
        }

        // now, column by column, add in the s*n contribution and reduce to 
        // a single 64+x bit accumulator while storing the intermediate product
        // 's' as we go.
        

#if 1
        j = 0;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te0, to0);
        acc_e1 = _mm512_add_epi64(acc_e1, te1);
        acc_o1 = _mm512_add_epi64(acc_o1, to1);

        prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
        prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
        a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

        // add in the final product
        ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 1;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te2, to2);
        acc_e1 = _mm512_add_epi64(acc_e1, te3);
        acc_o1 = _mm512_add_epi64(acc_o1, to3);

        for (k = 0; k < j; k++)
        {
            a0 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + k) * VECLEN);
            b0 = _mm512_load_epi32(n->data + (j - k) * VECLEN);

            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);
        }

        prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
        prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
        a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

        // add in the final product
        ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 2;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te4, to4);
        acc_e1 = _mm512_add_epi64(acc_e1, te5);
        acc_o1 = _mm512_add_epi64(acc_o1, to5);

        for (k = 0; k < j; k++)
        {
            a0 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + k) * VECLEN);
            b0 = _mm512_load_epi32(n->data + (j - k) * VECLEN);

            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);
        }

        prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
        prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
        a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

        // add in the final product
        ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 3;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te6, to6);
        acc_e1 = _mm512_add_epi64(acc_e1, te7);
        acc_o1 = _mm512_add_epi64(acc_o1, to7);

        for (k = 0; k < j; k++)
        {
            a0 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + k) * VECLEN);
            b0 = _mm512_load_epi32(n->data + (j - k) * VECLEN);

            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);
        }

        prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
        prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
        a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

        // add in the final product
        ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

#else

        for (j = 0; j < BLOCKWORDS; j++)
        {
            switch (j)
            {
            case 0:
                // accumulate this column-sum
                ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te0, to0);
                acc_e1 = _mm512_add_epi64(acc_e1, te1);
                acc_o1 = _mm512_add_epi64(acc_o1, to1);
                break;

            case 1:
                // accumulate this column-sum
                ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te2, to2);
                acc_e1 = _mm512_add_epi64(acc_e1, te3);
                acc_o1 = _mm512_add_epi64(acc_o1, to3);
                break;

            case 2:
                // accumulate this column-sum
                ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te4, to4);
                acc_e1 = _mm512_add_epi64(acc_e1, te5);
                acc_o1 = _mm512_add_epi64(acc_o1, to5);
                break;

            case 3:
                // accumulate this column-sum
                ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te6, to6);
                acc_e1 = _mm512_add_epi64(acc_e1, te7);
                acc_o1 = _mm512_add_epi64(acc_o1, to7);
                break;
            }

            for (k = 0; k < j; k++)
            {
                a0 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + k) * VECLEN);
                b0 = _mm512_load_epi32(n->data + (j - k) * VECLEN);

                _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
                ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);
            }

            prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
            prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
            a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

            _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

            b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
            _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

            // add in the final product
            ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

            // now shift.
            acc_e0 = _mm512_srli_epi64(acc_e0, 32);
            acc_o0 = _mm512_srli_epi64(acc_o0, 32);
            acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
            acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
            acc_e1 = zero;
            acc_o1 = zero;
        }
#endif
    }

    // second half mul
    for (i = NBLOCKS; i < 2 * NBLOCKS; i++)
    {

        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;
        to0 = to1 = to2 = to3 = to4 = to5 = to6 = to7 = zero;

        for (j = i - NBLOCKS + 1; j < NBLOCKS; j++)
        {
            a0 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);

            __m512i prod1_e, prod2_e, prod3_e, prod4_e, prod5_e;
            __mmask8 scarry_e1, scarry_e2, scarry_e3, scarry_e4;

            // a0 terms
            MUL_ACCUM_4X(a0, b0, b1, b2, b3, hiword);

            // a1 terms
            b0 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            MUL_ACCUM_4X(a1, b1, b2, b3, b0, hiword);

            // a2 terms
            b1 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            MUL_ACCUM_4X(a2, b2, b3, b0, b1, hiword);

            // a3 terms
            b2 = _mm512_load_epi32(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);
            MUL_ACCUM_4X(a3, b3, b0, b1, b2, hiword);
        }

        for (j = i - NBLOCKS + 1; j < NBLOCKS; j++)
        {
            // accumulate s * n
            a0 = _mm512_load_epi32(s->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi32(s->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi32(s->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi32(s->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);

            __m512i prod1_e, prod2_e, prod3_e, prod4_e, prod5_e;
            __mmask8 scarry_e1, scarry_e2, scarry_e3, scarry_e4;

            // a0 terms
            MUL_ACCUM_4X(a0, b0, b1, b2, b3, hiword);

            // a1 terms
            b0 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            MUL_ACCUM_4X(a1, b1, b2, b3, b0, hiword);

            // a2 terms
            b1 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            MUL_ACCUM_4X(a2, b2, b3, b0, b1, hiword);

            // a3 terms
            b2 = _mm512_load_epi32(n->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);
            MUL_ACCUM_4X(a3, b3, b0, b1, b2, hiword);
        }

        // finish each triangular shaped column sum (a * b)
        a1 = _mm512_load_epi32(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi32(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi32(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi32(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi32(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi32(b->data + (NWORDS - 3) * VECLEN);

        // ======
        _mm512_mul_eo64_epi32(a1, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        _mm512_mul_eo64_epi32(a2, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        _mm512_mul_eo64_epi32(a3, b2, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        // ======
        _mm512_mul_eo64_epi32(a2, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te2, to2, te3, to3);

        _mm512_mul_eo64_epi32(a3, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te2, to2, te3, to3);

        // ======
        _mm512_mul_eo64_epi32(a3, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te4, to4, te5, to5);

        // finish each triangular shaped column sum (s * n)
        a1 = _mm512_load_epi32(s->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi32(s->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi32(s->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi32(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi32(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi32(n->data + (NWORDS - 3) * VECLEN);

        // ======
        _mm512_mul_eo64_epi32(a1, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        _mm512_mul_eo64_epi32(a2, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        _mm512_mul_eo64_epi32(a3, b2, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        // ======
        _mm512_mul_eo64_epi32(a2, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te2, to2, te3, to3);

        _mm512_mul_eo64_epi32(a3, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te2, to2, te3, to3);

        // ======
        _mm512_mul_eo64_epi32(a3, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te4, to4, te5, to5);

        

#if 1
        j = 0;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te0, to0);
        acc_e1 = _mm512_add_epi64(acc_e1, te1);
        acc_o1 = _mm512_add_epi64(acc_o1, to1);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + ((i - NBLOCKS) * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 1;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te2, to2);
        acc_e1 = _mm512_add_epi64(acc_e1, te3);
        acc_o1 = _mm512_add_epi64(acc_o1, to3);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + ((i - NBLOCKS) * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 2;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te4, to4);
        acc_e1 = _mm512_add_epi64(acc_e1, te5);
        acc_o1 = _mm512_add_epi64(acc_o1, to5);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + ((i - NBLOCKS) * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 3;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te6, to6);
        acc_e1 = _mm512_add_epi64(acc_e1, te7);
        acc_o1 = _mm512_add_epi64(acc_o1, to7);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + ((i - NBLOCKS) * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

#else
        for (j = 0; j < BLOCKWORDS; j++)
        {
            switch (j)
            {
            case 0:
                // accumulate this column-sum
                ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te0, to0);
                acc_e1 = _mm512_add_epi64(acc_e1, te1);
                acc_o1 = _mm512_add_epi64(acc_o1, to1);
                break;

            case 1:
                // accumulate this column-sum
                ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te2, to2);
                acc_e1 = _mm512_add_epi64(acc_e1, te3);
                acc_o1 = _mm512_add_epi64(acc_o1, to3);
                break;

            case 2:
                // accumulate this column-sum
                ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te4, to4);
                acc_e1 = _mm512_add_epi64(acc_e1, te5);
                acc_o1 = _mm512_add_epi64(acc_o1, to5);
                break;

            case 3:
                // accumulate this column-sum
                ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te6, to6);
                acc_e1 = _mm512_add_epi64(acc_e1, te7);
                acc_o1 = _mm512_add_epi64(acc_o1, to7);
                break;
            }

            // store the low-word final result
            a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
            _mm512_store_epi32(s->data + ((i - NBLOCKS) * BLOCKWORDS + j) * VECLEN, a0);

            // and shift.
            acc_e0 = _mm512_srli_epi64(acc_e0, 32);
            acc_o0 = _mm512_srli_epi64(acc_o0, 32);
            acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
            acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
            acc_e1 = zero;
            acc_o1 = zero;
        }

#endif
    }

    a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
    
#if 0
    scarry2 = _mm512_cmp_epu32_mask(a0, zero, _MM_CMPINT_EQ);

    // subtract n from tmp
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a1 = _mm512_load_epi32(s->data + i * VECLEN);
        b0 = _mm512_load_epi32(n->data + i * VECLEN);
        a0 = _mm512_sbb_epi32(a1, scarry, b0, &scarry);
        _mm512_store_epi32(c->data + i * VECLEN, a0);
    }

    // negate any final borrows if there was also a final carry.
    scarry &= scarry2;

    // if there was a final borrow, we didn't need to do the subtraction after all.
    // replace with original results based on final borrow mask.
    for (i = NWORDS - 1; i >= 0; i--)
    {
        b0 = _mm512_load_epi32(s->data + i * VECLEN);
        _mm512_mask_store_epi32(c->data + i * VECLEN, scarry, b0);
    }

#else

    scarry2 = _mm512_cmp_epu32_mask(a0, zero, _MM_CMPINT_GT);
    scarry2 |= vec_mask_gte(scarry2, s, n);

    // when necessary, subtract n from the result
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a1 = _mm512_load_epi32(s->data + i * VECLEN);
        b0 = _mm512_load_epi32(n->data + i * VECLEN);
        a0 = _mm512_mask_sbb_epi32(a1, scarry2, scarry, b0, &scarry);
        _mm512_store_epi32(c->data + i * VECLEN, a0);
    }
    _mm512_store_epi32(c->data + i * VECLEN, _mm512_set1_epi32(0));

#endif

    c->size = NWORDS;
    return;
}

void vecsqrmod_bps(bignum* a, bignum* c, bignum* n, bignum* s, monty* mdata)
{
    int i, j, k;
    bignum* b = a;

    __m512i a0, a1, a2, a3;
    __m512i b0, b1, b2, b3, b4, b5, b6;
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;
    __m512i to0, to1, to2, to3, to4, to5, to6, to7;

    __m512i acc_e0;
    __m512i acc_o0;

    __m512i acc_e1;
    __m512i acc_o1;
    // 31

    __m512i nhatvec_e = _mm512_load_epi32(mdata->vrho); // _mm512_set1_epi32(nhat);
    __m512i nhatvec_o = _mm512_shuffle_epi32(nhatvec_e, 0xB1);;

    __m512i prod1_e;
    __m512i prod1_o;

    __m512i hiword = _mm512_set1_epi64(0x000000100000000);
    __m512i hiword2 = _mm512_set1_epi64(0x000000200000000);
    __m512i zero = _mm512_set1_epi64(0);
    // 37

    __mmask8 scarry_e1 = 0;
    __mmask8 scarry_o1 = 0;
    __mmask16 scarry2;
    __mmask16 scarry;

    // zero the accumulator
    acc_e0 = acc_o0 = acc_e1 = acc_o1 = zero;


    // first half sqr
    for (i = 0; i < NBLOCKS; i++)
    {

        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;
        to0 = to1 = to2 = to3 = to4 = to5 = to6 = to7 = zero;

        if (i & 1)
        {
            __mmask8 scarry_e1 = 0;
            __mmask8 scarry_o1 = 0;

            // i odd
            for (j = 0; j < (i - 1) / 2; j++)
            {
                a0 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 1) * VECLEN);
                a1 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 2) * VECLEN);
                a2 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 3) * VECLEN);
                a3 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 4) * VECLEN);

                b0 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 1) * VECLEN);
                b1 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 2) * VECLEN);
                b2 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 3) * VECLEN);
                b3 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 4) * VECLEN);

                __m512i prod1_e, prod2_e, prod3_e, prod4_e, prod5_e;
                __mmask8 scarry_e1, scarry_e2, scarry_e3, scarry_e4;

                // a0 terms
                MUL_ACCUM_4X(a0, b0, b1, b2, b3, hiword2);

                // a1 terms
                b0 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 5) * VECLEN);
                MUL_ACCUM_4X(a1, b1, b2, b3, b0, hiword2);

                // a2 terms
                b1 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 6) * VECLEN);
                MUL_ACCUM_4X(a2, b2, b3, b0, b1, hiword2);

                // a3 terms
                b2 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 7) * VECLEN);
                MUL_ACCUM_4X(a3, b3, b0, b1, b2, hiword2);
            }

            // i,j = (1,0), (3,1), (5,2), ... 
            a0 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS - 1) * VECLEN);  // 3,7,11...
            a1 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS - 2) * VECLEN);  // 2,6,10...
            a2 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS - 3) * VECLEN);  // 1,5,9...
            a3 = _mm512_load_epi32(a->data + ((i - j) * BLOCKWORDS - 4) * VECLEN);  // 0,4,8...
                                                                                    // for i=0,1,2...
            //b2 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + 3) * VECLEN);        // 3,7,11 (always = a0)
            b3 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + 4) * VECLEN);        // 4,8,12
            b4 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + 5) * VECLEN);        // 5,9,13
            b5 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + 6) * VECLEN);        // 6,10,14
            b6 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + 7) * VECLEN);        // 7,11,15

            //k == 0;
            _mm512_mul_eo64_epi32(a2, a0, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a3, b3, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

            //k == 1;
            _mm512_mul_eo64_epi32(a1, a0, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a3, b4, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            //k == 2;
            _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a2, b4, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a3, b5, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            //k == 3;
            _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a2, b5, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a3, b6, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            // the square (non-doubled) terms
            b2 = a0;
            a0 = a1;
            a1 = b2;
        }
        else
        {
            __mmask8 scarry_e1 = 0;
            __mmask8 scarry_o1 = 0;

            // i even
            for (j = 0; j < i / 2; j++)
            {
                a0 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 1) * VECLEN);
                a1 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 2) * VECLEN);
                a2 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 3) * VECLEN);
                a3 = _mm512_load_epi32(a->data + ((j + 1) * BLOCKWORDS - 4) * VECLEN);

                b0 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j
                    * BLOCKWORDS + 1) * VECLEN);
                b1 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j
                    * BLOCKWORDS + 2) * VECLEN);
                b2 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j
                    * BLOCKWORDS + 3) * VECLEN);
                b3 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j
                    * BLOCKWORDS + 4) * VECLEN);

                __m512i prod1_e, prod2_e, prod3_e, prod4_e, prod5_e;
                __mmask8 scarry_e1, scarry_e2, scarry_e3, scarry_e4;

                // a0 terms
                MUL_ACCUM_4X(a0, b0, b1, b2, b3, hiword2);

                // a1 terms
                b0 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 5) * VECLEN);
                MUL_ACCUM_4X(a1, b1, b2, b3, b0, hiword2);

                // a2 terms
                b1 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 6) * VECLEN);
                MUL_ACCUM_4X(a2, b2, b3, b0, b1, hiword2);

                // a3 terms
                b2 = _mm512_load_epi32(b->data + ((i - 1) * BLOCKWORDS - j * BLOCKWORDS + 7) * VECLEN);
                MUL_ACCUM_4X(a3, b3, b0, b1, b2, hiword2);
            }

            a0 = _mm512_load_epi32(a->data + (i / 2 * BLOCKWORDS + 0) * VECLEN);
            a1 = _mm512_load_epi32(a->data + (i / 2 * BLOCKWORDS + 1) * VECLEN);
            a2 = _mm512_load_epi32(a->data + (i / 2 * BLOCKWORDS + 2) * VECLEN);
            a3 = _mm512_load_epi32(a->data + (i / 2 * BLOCKWORDS + 3) * VECLEN);

            // save independent sum/carry words for each product-column in the block.
            // uses 11 register inputs, 16 register outputs, and 3 aux vectors.

            //k == 1;
            _mm512_mul_eo64_epi32(a0, a1, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            //k == 2;
            _mm512_mul_eo64_epi32(a0, a2, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            //k == 3;
            _mm512_mul_eo64_epi32(a0, a3, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a1, a2, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);
        }

        // all terms so far need to be doubled.  Do that all at once with these
        // left shifts.
        te1 = _mm512_or_epi64(te1, _mm512_maskz_srli_epi32(0xaaaa, te0, 31));
        to1 = _mm512_or_epi64(to1, _mm512_maskz_srli_epi32(0xaaaa, to0, 31));
        te0 = _mm512_slli_epi64(te0, 1);
        to0 = _mm512_slli_epi64(to0, 1);

        te3 = _mm512_or_epi64(te3, _mm512_maskz_srli_epi32(0xaaaa, te2, 31));
        to3 = _mm512_or_epi64(to3, _mm512_maskz_srli_epi32(0xaaaa, to2, 31));
        te2 = _mm512_slli_epi64(te2, 1);
        to2 = _mm512_slli_epi64(to2, 1);

        te5 = _mm512_or_epi64(te5, _mm512_maskz_srli_epi32(0xaaaa, te4, 31));
        to5 = _mm512_or_epi64(to5, _mm512_maskz_srli_epi32(0xaaaa, to4, 31));
        te4 = _mm512_slli_epi64(te4, 1);
        to4 = _mm512_slli_epi64(to4, 1);

        te7 = _mm512_or_epi64(te7, _mm512_maskz_srli_epi32(0xaaaa, te6, 31));
        to7 = _mm512_or_epi64(to7, _mm512_maskz_srli_epi32(0xaaaa, to6, 31));
        te6 = _mm512_slli_epi64(te6, 1);
        to6 = _mm512_slli_epi64(to6, 1);

        // finally the two non-doubled terms.
        _mm512_mul_eo64_epi32(a0, a0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        _mm512_mul_eo64_epi32(a1, a1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te4, to4, te5, to5);


        // for those 's' we have already accumulated, compute the
        // block s*n accumulations
        for (j = 0; j < i; j++)
        {
            // accumulate s * n
            a0 = _mm512_load_epi32(s->data + ((j + 1) * BLOCKWORDS - 1) * VECLEN);
            a1 = _mm512_load_epi32(s->data + ((j + 1) * BLOCKWORDS - 2) * VECLEN);
            a2 = _mm512_load_epi32(s->data + ((j + 1) * BLOCKWORDS - 3) * VECLEN);
            a3 = _mm512_load_epi32(s->data + ((j + 1) * BLOCKWORDS - 4) * VECLEN);

            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 4) * VECLEN);

            __m512i prod1_e, prod2_e, prod3_e, prod4_e, prod5_e;
            __mmask8 scarry_e1, scarry_e2, scarry_e3, scarry_e4;

            // a0 terms
            MUL_ACCUM_4X(a0, b0, b1, b2, b3, hiword);

            // a1 terms
            b0 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 5) * VECLEN);
            MUL_ACCUM_4X(a1, b1, b2, b3, b0, hiword);

            // a2 terms
            b1 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 6) * VECLEN);
            MUL_ACCUM_4X(a2, b2, b3, b0, b1, hiword);

            // a3 terms
            b2 = _mm512_load_epi32(n->data + ((i - j - 1) * BLOCKWORDS + 7) * VECLEN);
            MUL_ACCUM_4X(a3, b3, b0, b1, b2, hiword);
        }

#if 1
        j = 0;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te0, to0);
        acc_e1 = _mm512_add_epi64(acc_e1, te1);
        acc_o1 = _mm512_add_epi64(acc_o1, to1);

        prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
        prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
        a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

        // add in the final product
        ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 1;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te2, to2);
        acc_e1 = _mm512_add_epi64(acc_e1, te3);
        acc_o1 = _mm512_add_epi64(acc_o1, to3);

        for (k = 0; k < j; k++)
        {
            a0 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + k) * VECLEN);
            b0 = _mm512_load_epi32(n->data + (j - k) * VECLEN);

            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);
        }

        prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
        prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
        a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

        // add in the final product
        ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 2;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te4, to4);
        acc_e1 = _mm512_add_epi64(acc_e1, te5);
        acc_o1 = _mm512_add_epi64(acc_o1, to5);

        for (k = 0; k < j; k++)
        {
            a0 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + k) * VECLEN);
            b0 = _mm512_load_epi32(n->data + (j - k) * VECLEN);

            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);
        }

        prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
        prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
        a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

        // add in the final product
        ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 3;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te6, to6);
        acc_e1 = _mm512_add_epi64(acc_e1, te7);
        acc_o1 = _mm512_add_epi64(acc_o1, to7);

        for (k = 0; k < j; k++)
        {
            a0 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + k) * VECLEN);
            b0 = _mm512_load_epi32(n->data + (j - k) * VECLEN);

            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);
        }

        prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
        prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
        a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

        // add in the final product
        ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

#else

        for (j = 0; j < BLOCKWORDS; j++)
        {
            switch (j)
            {
            case 0:
                // accumulate this column-sum
                ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te0, to0);
                acc_e1 = _mm512_add_epi64(acc_e1, te1);
                acc_o1 = _mm512_add_epi64(acc_o1, to1);
                break;

            case 1:
                // accumulate this column-sum
                ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te2, to2);
                acc_e1 = _mm512_add_epi64(acc_e1, te3);
                acc_o1 = _mm512_add_epi64(acc_o1, to3);
                break;

            case 2:
                // accumulate this column-sum
                ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te4, to4);
                acc_e1 = _mm512_add_epi64(acc_e1, te5);
                acc_o1 = _mm512_add_epi64(acc_o1, to5);
                break;

            case 3:
                // accumulate this column-sum
                ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te6, to6);
                acc_e1 = _mm512_add_epi64(acc_e1, te7);
                acc_o1 = _mm512_add_epi64(acc_o1, to7);
                break;
            }

            for (k = 0; k < j; k++)
            {
                a0 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + k) * VECLEN);
                b0 = _mm512_load_epi32(n->data + (j - k) * VECLEN);

                _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
                ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);
            }

            prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
            prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
            a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

            _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

            b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
            _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

            // add in the final product
            ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

            // now shift.
            acc_e0 = _mm512_srli_epi64(acc_e0, 32);
            acc_o0 = _mm512_srli_epi64(acc_o0, 32);
            acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
            acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
            acc_e1 = zero;
            acc_o1 = zero;
        }
#endif
    }

    // second half sqr
    for (i = 0; i < NBLOCKS; i++)
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;
        to0 = to1 = to2 = to3 = to4 = to5 = to6 = to7 = zero;

        for (j = 0; j < (NBLOCKS - i - 1) / 2; j++)
        {
            // Compute a solid block (all matching terms are in the lower
            // half triangle of the expansion).

            //hips_block_mul_type3(a->data + (i * BLOCKWORDS + (j + 2) * BLOCKWORDS) * VECLEN,
            //    a->data + (NWORDS - 2 * BLOCKWORDS - j * BLOCKWORDS) * VECLEN, t_e, t_o);
            a0 = _mm512_load_epi32(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
            a1 = _mm512_load_epi32(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
            a2 = _mm512_load_epi32(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);
            a3 = _mm512_load_epi32(a->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);

            b0 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 4) * VECLEN);

            __m512i prod1_e, prod2_e, prod3_e, prod4_e, prod5_e;
            __mmask8 scarry_e1, scarry_e2, scarry_e3, scarry_e4;

            // a0 terms
            MUL_ACCUM_4X(a0, b0, b1, b2, b3, hiword2);

            // a1 terms
            b0 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 5) * VECLEN);
            MUL_ACCUM_4X(a1, b1, b2, b3, b0, hiword2);

            // a2 terms
            b1 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 6) * VECLEN);
            MUL_ACCUM_4X(a2, b2, b3, b0, b1, hiword2);

            // a3 terms
            b2 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 7) * VECLEN);
            MUL_ACCUM_4X(a3, b3, b0, b1, b2, hiword2);
        }

        // The final block shape depends on the parity of i and NBLOCKS
        // Block shape 1 (small upper triangle) if 'i' is odd and NBLOCKS is even,
        // or if 'i' is even and NBLOCKS is odd.
        // Block shape 2 if 'i' is odd and NBLOCKS is odd or if 'i' is even
        // and NBLOCKS is even.
#if NBLOCKS & 1		// NBLOCKS is odd
// i odd, block shape 2.
        if (i & 1)
        {
            // always a continuation of the full-block loop, so use the same 
            // loading pattern.  Only now we don't need as many b-terms.
            a0 = _mm512_load_epi32(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
            a1 = _mm512_load_epi32(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
            a2 = _mm512_load_epi32(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);
            a3 = _mm512_load_epi32(a->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);

            b0 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 5) * VECLEN);

            // save independent sum/carry words for each product-column in the block.
            // uses 11 register inputs, 16 register outputs, and 3 aux vectors.
            //k == 0;
            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

            //k == 1;
            _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            //k == 2;
            _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            //k == 3;
            _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            // all terms so far need to be doubled.  Do that all at once with these
            // left shifts.
            te1 = _mm512_or_epi64(te1, _mm512_maskz_srli_epi32(0xaaaa, te0, 31));
            to1 = _mm512_or_epi64(to1, _mm512_maskz_srli_epi32(0xaaaa, to0, 31));
            te0 = _mm512_slli_epi64(te0, 1);
            to0 = _mm512_slli_epi64(to0, 1);

            te3 = _mm512_or_epi64(te3, _mm512_maskz_srli_epi32(0xaaaa, te2, 31));
            to3 = _mm512_or_epi64(to3, _mm512_maskz_srli_epi32(0xaaaa, to2, 31));
            te2 = _mm512_slli_epi64(te2, 1);
            to2 = _mm512_slli_epi64(to2, 1);

            te5 = _mm512_or_epi64(te5, _mm512_maskz_srli_epi32(0xaaaa, te4, 31));
            to5 = _mm512_or_epi64(to5, _mm512_maskz_srli_epi32(0xaaaa, to4, 31));
            te4 = _mm512_slli_epi64(te4, 1);
            to4 = _mm512_slli_epi64(to4, 1);

            te7 = _mm512_or_epi64(te7, _mm512_maskz_srli_epi32(0xaaaa, te6, 31));
            to7 = _mm512_or_epi64(to7, _mm512_maskz_srli_epi32(0xaaaa, to6, 31));
            te6 = _mm512_slli_epi64(te6, 1);
            to6 = _mm512_slli_epi64(to6, 1);

            // finally the two non-doubled terms.
            _mm512_mul_eo64_epi32(a3, a3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a2, a2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);
        }
        else
        {
            // i even, block shape 1.
            // always a continuation of the full-block loop, so use the same 
            // loading pattern.  Only now we don't need as many b-terms.
            a0 = _mm512_load_epi32(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
            a1 = _mm512_load_epi32(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
            a2 = _mm512_load_epi32(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);

            //k == 0;
            _mm512_mul_eo64_epi32(a0, a2, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

            //k == 1;
            _mm512_mul_eo64_epi32(a0, a1, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            te1 = _mm512_or_epi64(te1, _mm512_maskz_srli_epi32(0xaaaa, te0, 31));
            to1 = _mm512_or_epi64(to1, _mm512_maskz_srli_epi32(0xaaaa, to0, 31));
            te0 = _mm512_slli_epi64(te0, 1);
            to0 = _mm512_slli_epi64(to0, 1);

            te3 = _mm512_or_epi64(te3, _mm512_maskz_srli_epi32(0xaaaa, te2, 31));
            to3 = _mm512_or_epi64(to3, _mm512_maskz_srli_epi32(0xaaaa, to2, 31));
            te2 = _mm512_slli_epi64(te2, 1);
            to2 = _mm512_slli_epi64(to2, 1);

            // technically only have to do these two if j > 0 (so that
            // they are non-zero from full-block loop iterations).
            te5 = _mm512_or_epi64(te5, _mm512_maskz_srli_epi32(0xaaaa, te4, 31));
            to5 = _mm512_or_epi64(to5, _mm512_maskz_srli_epi32(0xaaaa, to4, 31));
            te4 = _mm512_slli_epi64(te4, 1);
            to4 = _mm512_slli_epi64(to4, 1);

            te7 = _mm512_or_epi64(te7, _mm512_maskz_srli_epi32(0xaaaa, te6, 31));
            to7 = _mm512_or_epi64(to7, _mm512_maskz_srli_epi32(0xaaaa, to6, 31));
            te6 = _mm512_slli_epi64(te6, 1);
            to6 = _mm512_slli_epi64(to6, 1);

            // finally the two non-doubled terms.
            _mm512_mul_eo64_epi32(a1, a1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a0, a0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);
        }


#else				// NBLOCKS is even
// i odd, block shape 1.
        if (i & 1)
        {
            // always a continuation of the full-block loop, so use the same 
            // loading pattern.  Only now we don't need as many b-terms.
            a0 = _mm512_load_epi32(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);  // {f, b}
            a1 = _mm512_load_epi32(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);  // {e, a}
            a2 = _mm512_load_epi32(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);  // {d, 9}

            //k == 0;
            _mm512_mul_eo64_epi32(a0, a2, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

            //k == 1;
            _mm512_mul_eo64_epi32(a0, a1, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            te1 = _mm512_or_epi64(te1, _mm512_maskz_srli_epi32(0xaaaa, te0, 31));
            to1 = _mm512_or_epi64(to1, _mm512_maskz_srli_epi32(0xaaaa, to0, 31));
            te0 = _mm512_slli_epi64(te0, 1);
            to0 = _mm512_slli_epi64(to0, 1);

            te3 = _mm512_or_epi64(te3, _mm512_maskz_srli_epi32(0xaaaa, te2, 31));
            to3 = _mm512_or_epi64(to3, _mm512_maskz_srli_epi32(0xaaaa, to2, 31));
            te2 = _mm512_slli_epi64(te2, 1);
            to2 = _mm512_slli_epi64(to2, 1);

            // technically only have to do these two if j > 0 (so that
            // they are non-zero from full-block loop iterations).
            te5 = _mm512_or_epi64(te5, _mm512_maskz_srli_epi32(0xaaaa, te4, 31));
            to5 = _mm512_or_epi64(to5, _mm512_maskz_srli_epi32(0xaaaa, to4, 31));
            te4 = _mm512_slli_epi64(te4, 1);
            to4 = _mm512_slli_epi64(to4, 1);

            te7 = _mm512_or_epi64(te7, _mm512_maskz_srli_epi32(0xaaaa, te6, 31));
            to7 = _mm512_or_epi64(to7, _mm512_maskz_srli_epi32(0xaaaa, to6, 31));
            te6 = _mm512_slli_epi64(te6, 1);
            to6 = _mm512_slli_epi64(to6, 1);

            // finally the two non-doubled terms.
            _mm512_mul_eo64_epi32(a1, a1, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a0, a0, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);
        }
        else
        {
            // i even, block shape 1.
            // always a continuation of the full-block loop, so use the same 
            // loading pattern.  Only now we don't need as many b-terms.
            a0 = _mm512_load_epi32(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);		// {f, b}
            a1 = _mm512_load_epi32(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);		// {e, a}
            a2 = _mm512_load_epi32(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);		// {d, 9}
            a3 = _mm512_load_epi32(a->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);		// {c, 8}

            b0 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 1) * VECLEN); // {9, 5}
            b1 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 2) * VECLEN);	// {a, 6}
            b2 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 3) * VECLEN);	// {b, 7}
            b3 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 4) * VECLEN);	// {c, 8}
            b4 = _mm512_load_epi32(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 5) * VECLEN); // {d, 9}

            // save independent sum/carry words for each product-column in the block.
            // uses 11 register inputs, 16 register outputs, and 3 aux vectors.
            //k == 0;
            _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a1, b1, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a2, b2, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te0, to0, te1, to1);

            //k == 1;
            _mm512_mul_eo64_epi32(a0, b1, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a1, b2, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            _mm512_mul_eo64_epi32(a2, b3, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te2, to2, te3, to3);

            //k == 2;
            _mm512_mul_eo64_epi32(a0, b2, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            _mm512_mul_eo64_epi32(a1, b3, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te4, to4, te5, to5);

            //k == 3;
            _mm512_mul_eo64_epi32(a0, b3, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            _mm512_mul_eo64_epi32(a1, b4, &prod1_e, &prod1_o);
            ACCUM_DOUBLED_EO_PROD(te6, to6, te7, to7);

            // all terms so far need to be doubled.  Do that all at once with these
            // left shifts.
            te1 = _mm512_or_epi64(te1, _mm512_maskz_srli_epi32(0xaaaa, te0, 31));
            to1 = _mm512_or_epi64(to1, _mm512_maskz_srli_epi32(0xaaaa, to0, 31));
            te0 = _mm512_slli_epi64(te0, 1);
            to0 = _mm512_slli_epi64(to0, 1);

            te3 = _mm512_or_epi64(te3, _mm512_maskz_srli_epi32(0xaaaa, te2, 31));
            to3 = _mm512_or_epi64(to3, _mm512_maskz_srli_epi32(0xaaaa, to2, 31));
            te2 = _mm512_slli_epi64(te2, 1);
            to2 = _mm512_slli_epi64(to2, 1);

            te5 = _mm512_or_epi64(te5, _mm512_maskz_srli_epi32(0xaaaa, te4, 31));
            to5 = _mm512_or_epi64(to5, _mm512_maskz_srli_epi32(0xaaaa, to4, 31));
            te4 = _mm512_slli_epi64(te4, 1);
            to4 = _mm512_slli_epi64(to4, 1);

            te7 = _mm512_or_epi64(te7, _mm512_maskz_srli_epi32(0xaaaa, te6, 31));
            to7 = _mm512_or_epi64(to7, _mm512_maskz_srli_epi32(0xaaaa, to6, 31));
            te6 = _mm512_slli_epi64(te6, 1);
            to6 = _mm512_slli_epi64(to6, 1);

            // finally the two non-doubled terms.
            _mm512_mul_eo64_epi32(a3, a3, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te0, to0, te1, to1);

            _mm512_mul_eo64_epi32(a2, a2, &prod1_e, &prod1_o);
            ACCUM_EO_PROD(te4, to4, te5, to5);
        }
#endif

        // the s*n term.  No more doubling past here.
        for (j = 0; j < NBLOCKS - 1 - i; j++)
        {
            a0 = _mm512_load_epi32(s->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
            a1 = _mm512_load_epi32(s->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
            a2 = _mm512_load_epi32(s->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);
            a3 = _mm512_load_epi32(s->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);

            b0 = _mm512_load_epi32(n->data + ((i + j) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi32(n->data + ((i + j) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi32(n->data + ((i + j) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi32(n->data + ((i + j) * BLOCKWORDS + 4) * VECLEN);

            __m512i prod1_e, prod2_e, prod3_e, prod4_e, prod5_e;
            __mmask8 scarry_e1, scarry_e2, scarry_e3, scarry_e4;

            // a0 terms
            MUL_ACCUM_4X(a0, b0, b1, b2, b3, hiword);

            // a1 terms
            b0 = _mm512_load_epi32(n->data + ((i + j) * BLOCKWORDS + 5) * VECLEN);
            MUL_ACCUM_4X(a1, b1, b2, b3, b0, hiword);

            // a2 terms
            b1 = _mm512_load_epi32(n->data + ((i + j) * BLOCKWORDS + 6) * VECLEN);
            MUL_ACCUM_4X(a2, b2, b3, b0, b1, hiword);

            // a3 terms
            b2 = _mm512_load_epi32(n->data + ((i + j) * BLOCKWORDS + 7) * VECLEN);
            MUL_ACCUM_4X(a3, b3, b0, b1, b2, hiword);
        }

        // finish each triangluar shaped column sum (s * n)
        a1 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi32(s->data + (i * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi32(n->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi32(n->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi32(n->data + (NWORDS - 3) * VECLEN);

        // ======
        _mm512_mul_eo64_epi32(a1, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        _mm512_mul_eo64_epi32(a2, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        _mm512_mul_eo64_epi32(a3, b2, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        // ======
        _mm512_mul_eo64_epi32(a2, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te2, to2, te3, to3);

        _mm512_mul_eo64_epi32(a3, b1, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te2, to2, te3, to3);

        // ======
        _mm512_mul_eo64_epi32(a3, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te4, to4, te5, to5);

#if 1
        j = 0;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te0, to0);
        acc_e1 = _mm512_add_epi64(acc_e1, te1);
        acc_o1 = _mm512_add_epi64(acc_o1, to1);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 1;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te2, to2);
        acc_e1 = _mm512_add_epi64(acc_e1, te3);
        acc_o1 = _mm512_add_epi64(acc_o1, to3);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 2;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te4, to4);
        acc_e1 = _mm512_add_epi64(acc_e1, te5);
        acc_o1 = _mm512_add_epi64(acc_o1, to5);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

        j = 3;
        // accumulate this column-sum
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te6, to6);
        acc_e1 = _mm512_add_epi64(acc_e1, te7);
        acc_o1 = _mm512_add_epi64(acc_o1, to7);

        // store the low-word final result
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

        // and shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;

#else
        for (j = 0; j < BLOCKWORDS; j++)
        {
            switch (j)
            {
            case 0:
                // accumulate this column-sum
                ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te0, to0);
                acc_e1 = _mm512_add_epi64(acc_e1, te1);
                acc_o1 = _mm512_add_epi64(acc_o1, to1);
                break;

            case 1:
                // accumulate this column-sum
                ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te2, to2);
                acc_e1 = _mm512_add_epi64(acc_e1, te3);
                acc_o1 = _mm512_add_epi64(acc_o1, to3);
                break;

            case 2:
                // accumulate this column-sum
                ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te4, to4);
                acc_e1 = _mm512_add_epi64(acc_e1, te5);
                acc_o1 = _mm512_add_epi64(acc_o1, to5);
                break;

            case 3:
                // accumulate this column-sum
                ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te6, to6);
                acc_e1 = _mm512_add_epi64(acc_e1, te7);
                acc_o1 = _mm512_add_epi64(acc_o1, to7);
                break;
            }

            // store the low-word final result
            a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
            _mm512_store_epi32(s->data + (i * BLOCKWORDS + j) * VECLEN, a0);

            // and shift.
            acc_e0 = _mm512_srli_epi64(acc_e0, 32);
            acc_o0 = _mm512_srli_epi64(acc_o0, 32);
            acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
            acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
            acc_e1 = zero;
            acc_o1 = zero;
        }

#endif
    }

    a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);

#if 0
    scarry2 = _mm512_cmp_epu32_mask(a0, zero, _MM_CMPINT_EQ);

    // subtract n from tmp
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a1 = _mm512_load_epi32(s->data + i * VECLEN);
        b0 = _mm512_load_epi32(n->data + i * VECLEN);
        a0 = _mm512_sbb_epi32(a1, scarry, b0, &scarry);
        _mm512_store_epi32(c->data + i * VECLEN, a0);
    }

    // negate any final borrows if there was also a final carry.
    scarry &= scarry2;

    // if there was a final borrow, we didn't need to do the subtraction after all.
    // replace with original results based on final borrow mask.
    for (i = NWORDS - 1; i >= 0; i--)
    {
        b0 = _mm512_load_epi32(s->data + i * VECLEN);
        _mm512_mask_store_epi32(c->data + i * VECLEN, scarry, b0);
    }

#else

    scarry2 = _mm512_cmp_epu32_mask(a0, zero, _MM_CMPINT_GT);
    scarry2 |= vec_mask_gte(scarry2, s, n);

    // when necessary, subtract n from the result
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a1 = _mm512_load_epi32(s->data + i * VECLEN);
        b0 = _mm512_load_epi32(n->data + i * VECLEN);
        a0 = _mm512_mask_sbb_epi32(a1, scarry2, scarry, b0, &scarry);
        _mm512_store_epi32(c->data + i * VECLEN, a0);
    }
    _mm512_store_epi32(c->data + i * VECLEN, _mm512_set1_epi32(0));

#endif

    c->size = NWORDS;
    return;
}
//_fips
void vecmulmod(bignum* a, bignum* b, bignum* c, bignum* n, bignum* s, monty* mdata)
{
    int i, j;

    __m512i a0, b0, s0, n0;
    __m512i te0, te1;
    __m512i to0, to1;

    __m512i acc_e0;
    __m512i acc_o0;

    __m512i acc_e1;
    __m512i acc_o1;

    __m512i nhatvec_e = _mm512_load_epi32(mdata->vrho);
    __m512i nhatvec_o = _mm512_shuffle_epi32(nhatvec_e, 0xB1);;

    __m512i prod1;
    __m512i prod2;

    __m512i prod1_e;
    __m512i prod1_o;

    __m512i hiword = _mm512_set1_epi64(0x000000100000000);
    __m512i zero = _mm512_set1_epi64(0);

    __mmask8 scarry_1 = 0;
    __mmask8 scarry_2 = 0;
    __mmask8 scarry_e1 = 0;
    __mmask8 scarry_o1 = 0;
    __mmask16 scarry2;
    __mmask16 scarry;

    // zero the accumulator
    acc_e0 = acc_o0 = acc_e1 = acc_o1 = zero;

    // first half mul
    for (i = 0; i < NWORDS; i++)
    {
        te0 = te1 = to0 = to1 = zero;
        for (j = 0; j < i; j++)
        {
            a0 = _mm512_load_epi32(a->data + j * VECLEN);
            b0 = _mm512_load_epi32(b->data + (i - j) * VECLEN);
            s0 = _mm512_load_epi32(s->data + j * VECLEN);
            n0 = _mm512_load_epi32(mdata->n->data + (i - j) * VECLEN);
            
            MUL_ACCUM_2X(a0, b0, s0, n0, hiword);
        }

        //spMulAddc(u->data[i], v->data[0], t);
        a0 = _mm512_load_epi32(a->data + i * VECLEN);
        b0 = _mm512_load_epi32(b->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(a0, b0, &prod1_e, &prod1_o);
        ACCUM_EO_PROD(te0, to0, te1, to1);

        // accumulate into main accumulator
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te0, to0);
        acc_e1 = _mm512_add_epi64(acc_e1, te1);
        acc_o1 = _mm512_add_epi64(acc_o1, to1);

        prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
        prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);
        a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);

        _mm512_store_epi32(s->data + i * VECLEN, a0);

        b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

        // add in the final product
        ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;
    }

    // second half mul
    for (i = NWORDS; i < 2 * NWORDS; i++)
    {
        te0 = te1 = to0 = to1 = zero;
        for (j = i - NWORDS + 1; j < NWORDS; j++)
        {
            a0 = _mm512_load_epi32(a->data + j * VECLEN);
            b0 = _mm512_load_epi32(b->data + (i - j) * VECLEN);
            s0 = _mm512_load_epi32(s->data + j * VECLEN);
            n0 = _mm512_load_epi32(mdata->n->data + (i - j) * VECLEN);

            MUL_ACCUM_2X(a0, b0, s0, n0, hiword);
            //spMul2Acc(u->data[j], v->data[i - j], mdata->n->data[i - j], s->data[j], t);
        }

        // accumulate into main accumulator
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te0, to0);
        acc_e1 = _mm512_add_epi64(acc_e1, te1);
        acc_o1 = _mm512_add_epi64(acc_o1, to1);

        // store
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + (i - NWORDS) * VECLEN, a0);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;
    }

    a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);

    scarry2 = _mm512_cmp_epu32_mask(a0, zero, _MM_CMPINT_GT);
    scarry2 |= vec_mask_gte(scarry2, s, n);

    // when necessary, subtract n from the result
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a0 = _mm512_load_epi32(s->data + i * VECLEN);
        b0 = _mm512_load_epi32(n->data + i * VECLEN);
        a0 = _mm512_mask_sbb_epi32(a0, scarry2, scarry, b0, &scarry);
        _mm512_store_epi32(c->data + i * VECLEN, a0);
    }
    _mm512_store_epi32(c->data + i * VECLEN, _mm512_set1_epi32(0));

    c->size = NWORDS;
    return;
}

void vecsqrmod(bignum* a, bignum* c, bignum* n, bignum* s, monty* mdata)
{
    int i, j;
    bignum* b = a;

    __m512i a0, b0;
    __m512i te0, te1;
    __m512i to0, to1;

    __m512i acc_e0;
    __m512i acc_o0;

    __m512i acc_e1;
    __m512i acc_o1;

    __m512i nhatvec_e = _mm512_load_epi32(mdata->vrho);
    __m512i nhatvec_o = _mm512_shuffle_epi32(nhatvec_e, 0xB1);;

    __m512i prod1;

    __m512i prod1_e;
    __m512i prod1_o;

    __m512i hiword = _mm512_set1_epi64(0x000000100000000);
    __m512i hiword2 = _mm512_set1_epi64(0x000000200000000);
    __m512i zero = _mm512_set1_epi64(0);

    __mmask8 scarry_1 = 0;
    __mmask8 scarry_e1 = 0;
    __mmask8 scarry_o1 = 0;
    __mmask16 scarry2;
    __mmask16 scarry;

    // zero the accumulator
    acc_e0 = acc_o0 = acc_e1 = acc_o1 = zero;

    // first half sqr
    for (i = 0; i < NWORDS; i++)
    {
        te0 = te1 = to0 = to1 = zero;

        // terms that are doubled
        for (j = 0; j < i / 2; j++)
        {
            a0 = _mm512_load_epi32(a->data + j * VECLEN);
            b0 = _mm512_load_epi32(b->data + (i - j) * VECLEN);

            prod1 = _mm512_mul_epu32(a0, b0);                                                                                
            te0 = _mm512_add_epi64(te0, prod1);                         
            scarry_1 = _mm512_cmplt_epu64_mask(te0, prod1);             
            te1 = _mm512_mask_add_epi64(te1, scarry_1, hiword2, te1);    
                                                                    
            a0 = _mm512_shuffle_epi32(a0, 0xB1);                    
            b0 = _mm512_shuffle_epi32(b0, 0xB1);                    
               
            prod1 = _mm512_mul_epu32(a0, b0);                       
            to0 = _mm512_add_epi64(to0, prod1);                     
            scarry_1 = _mm512_cmplt_epu64_mask(to0, prod1);         
            to1 = _mm512_mask_add_epi64(to1, scarry_1, hiword2, to1);
        }

        if ((i & 1) == 1)
        {
            // last doubled term
            a0 = _mm512_load_epi32(a->data + j * VECLEN);
            b0 = _mm512_load_epi32(b->data + (i - j) * VECLEN);

            prod1 = _mm512_mul_epu32(a0, b0);
            te0 = _mm512_add_epi64(te0, prod1);
            scarry_1 = _mm512_cmplt_epu64_mask(te0, prod1);
            te1 = _mm512_mask_add_epi64(te1, scarry_1, hiword2, te1);

            a0 = _mm512_shuffle_epi32(a0, 0xB1);
            b0 = _mm512_shuffle_epi32(b0, 0xB1);

            prod1 = _mm512_mul_epu32(a0, b0);
            to0 = _mm512_add_epi64(to0, prod1);
            scarry_1 = _mm512_cmplt_epu64_mask(to0, prod1);
            to1 = _mm512_mask_add_epi64(to1, scarry_1, hiword2, to1);
        }

        // double the column
        te1 = _mm512_or_epi64(te1, _mm512_maskz_srli_epi32(0xaaaa, te0, 31));
        to1 = _mm512_or_epi64(to1, _mm512_maskz_srli_epi32(0xaaaa, to0, 31));
        te0 = _mm512_slli_epi64(te0, 1);
        to0 = _mm512_slli_epi64(to0, 1);

        // s * n (not doubled)
        for (j = 0; j < i; j++)
        {
            a0 = _mm512_load_epi32(s->data + j * VECLEN);
            b0 = _mm512_load_epi32(mdata->n->data + (i - j) * VECLEN);

            prod1 = _mm512_mul_epu32(a0, b0);
            te0 = _mm512_add_epi64(te0, prod1);
            scarry_1 = _mm512_cmplt_epu64_mask(te0, prod1);
            te1 = _mm512_mask_add_epi64(te1, scarry_1, hiword, te1);

            a0 = _mm512_shuffle_epi32(a0, 0xB1);
            b0 = _mm512_shuffle_epi32(b0, 0xB1);

            prod1 = _mm512_mul_epu32(a0, b0);
            to0 = _mm512_add_epi64(to0, prod1);
            scarry_1 = _mm512_cmplt_epu64_mask(to0, prod1);
            to1 = _mm512_mask_add_epi64(to1, scarry_1, hiword, to1);
        }

        // square term (not doubled)
        if ((i & 1) == 0)
        {
            a0 = _mm512_load_epi32(a->data + i / 2 * VECLEN);

            prod1 = _mm512_mul_epu32(a0, a0);
            te0 = _mm512_add_epi64(te0, prod1);
            scarry_1 = _mm512_cmplt_epu64_mask(te0, prod1);
            te1 = _mm512_mask_add_epi64(te1, scarry_1, hiword, te1);

            a0 = _mm512_shuffle_epi32(a0, 0xB1);

            prod1 = _mm512_mul_epu32(a0, a0);
            to0 = _mm512_add_epi64(to0, prod1);
            scarry_1 = _mm512_cmplt_epu64_mask(to0, prod1);
            to1 = _mm512_mask_add_epi64(to1, scarry_1, hiword, to1);
        }

        // accumulate into main accumulator
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te0, to0);
        acc_e1 = _mm512_add_epi64(acc_e1, te1);
        acc_o1 = _mm512_add_epi64(acc_o1, to1);

        // mul with -n^-1 term
        prod1_e = _mm512_mul_epu32(nhatvec_e, acc_e0);
        prod1_o = _mm512_mul_epu32(nhatvec_o, acc_o0);

        // combine and store
        a0 = _mm512_eo64lo_to_epi32(prod1_e, prod1_o);
        _mm512_store_epi32(s->data + i * VECLEN, a0);

        // last s * n term
        b0 = _mm512_load_epi32(n->data + 0 * VECLEN);
        _mm512_mul_eo64_epi32(b0, a0, &prod1_e, &prod1_o);

        // add in the final product
        ACCUM_EO_PROD(acc_e0, acc_o0, acc_e1, acc_o1);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;
    }

    // second half mul
    for (i = NWORDS; i < 2 * NWORDS - 1; i++)
    {
        te0 = te1 = to0 = to1 = zero;
        
        // terms that are doubled
        for (j = i - NWORDS + 1; j < i / 2; j++)
        {
            a0 = _mm512_load_epi32(a->data + j * VECLEN);
            b0 = _mm512_load_epi32(b->data + (i - j) * VECLEN);

            prod1 = _mm512_mul_epu32(a0, b0);
            te0 = _mm512_add_epi64(te0, prod1);
            scarry_1 = _mm512_cmplt_epu64_mask(te0, prod1);
            te1 = _mm512_mask_add_epi64(te1, scarry_1, hiword2, te1);

            a0 = _mm512_shuffle_epi32(a0, 0xB1);
            b0 = _mm512_shuffle_epi32(b0, 0xB1);

            prod1 = _mm512_mul_epu32(a0, b0);
            to0 = _mm512_add_epi64(to0, prod1);
            scarry_1 = _mm512_cmplt_epu64_mask(to0, prod1);
            to1 = _mm512_mask_add_epi64(to1, scarry_1, hiword2, to1);
        }

        if (i & 1)
        {
            // last doubled term
            a0 = _mm512_load_epi32(a->data + j * VECLEN);
            b0 = _mm512_load_epi32(b->data + (i - j) * VECLEN);

            prod1 = _mm512_mul_epu32(a0, b0);
            te0 = _mm512_add_epi64(te0, prod1);
            scarry_1 = _mm512_cmplt_epu64_mask(te0, prod1);
            te1 = _mm512_mask_add_epi64(te1, scarry_1, hiword2, te1);

            a0 = _mm512_shuffle_epi32(a0, 0xB1);
            b0 = _mm512_shuffle_epi32(b0, 0xB1);

            prod1 = _mm512_mul_epu32(a0, b0);
            to0 = _mm512_add_epi64(to0, prod1);
            scarry_1 = _mm512_cmplt_epu64_mask(to0, prod1);
            to1 = _mm512_mask_add_epi64(to1, scarry_1, hiword2, to1);
        }

        // double the column
        te1 = _mm512_or_epi64(te1, _mm512_maskz_srli_epi32(0xaaaa, te0, 31));
        to1 = _mm512_or_epi64(to1, _mm512_maskz_srli_epi32(0xaaaa, to0, 31));
        te0 = _mm512_slli_epi64(te0, 1);
        to0 = _mm512_slli_epi64(to0, 1);

        // s * n (not doubled)
        for (j = i - NWORDS + 1; j < NWORDS; j++)
        {
            a0 = _mm512_load_epi32(s->data + j * VECLEN);
            b0 = _mm512_load_epi32(mdata->n->data + (i - j) * VECLEN);

            prod1 = _mm512_mul_epu32(a0, b0);
            te0 = _mm512_add_epi64(te0, prod1);
            scarry_1 = _mm512_cmplt_epu64_mask(te0, prod1);
            te1 = _mm512_mask_add_epi64(te1, scarry_1, hiword, te1);

            a0 = _mm512_shuffle_epi32(a0, 0xB1);
            b0 = _mm512_shuffle_epi32(b0, 0xB1);

            prod1 = _mm512_mul_epu32(a0, b0);
            to0 = _mm512_add_epi64(to0, prod1);
            scarry_1 = _mm512_cmplt_epu64_mask(to0, prod1);
            to1 = _mm512_mask_add_epi64(to1, scarry_1, hiword, to1);
        }

        // square term (not doubled)
        if ((i & 1) == 0)
        {
            a0 = _mm512_load_epi32(a->data + i / 2 * VECLEN);

            prod1 = _mm512_mul_epu32(a0, a0);
            te0 = _mm512_add_epi64(te0, prod1);
            scarry_1 = _mm512_cmplt_epu64_mask(te0, prod1);
            te1 = _mm512_mask_add_epi64(te1, scarry_1, hiword, te1);

            a0 = _mm512_shuffle_epi32(a0, 0xB1);

            prod1 = _mm512_mul_epu32(a0, a0);
            to0 = _mm512_add_epi64(to0, prod1);
            scarry_1 = _mm512_cmplt_epu64_mask(to0, prod1);
            to1 = _mm512_mask_add_epi64(to1, scarry_1, hiword, to1);
        }

        // accumulate into main accumulator
        ACCUM_EO_PROD2(acc_e0, acc_o0, acc_e1, acc_o1, te0, to0);
        acc_e1 = _mm512_add_epi64(acc_e1, te1);
        acc_o1 = _mm512_add_epi64(acc_o1, to1);

        // combine and store
        a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
        _mm512_store_epi32(s->data + (i - NWORDS) * VECLEN, a0);

        // now shift.
        acc_e0 = _mm512_srli_epi64(acc_e0, 32);
        acc_o0 = _mm512_srli_epi64(acc_o0, 32);
        acc_e0 = _mm512_add_epi64(acc_e1, acc_e0);
        acc_o0 = _mm512_add_epi64(acc_o1, acc_o0);
        acc_e1 = zero;
        acc_o1 = zero;
    }
    a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);
    _mm512_store_epi32(s->data + (i - NWORDS) * VECLEN, a0);

    acc_e0 = _mm512_srli_epi64(acc_e0, 32);
    acc_o0 = _mm512_srli_epi64(acc_o0, 32);
    a0 = _mm512_eo64lo_to_epi32(acc_e0, acc_o0);

    scarry2 = _mm512_cmp_epu32_mask(a0, zero, _MM_CMPINT_GT);
    scarry2 |= vec_mask_gte(scarry2, s, n);

    // when necessary, subtract n from the result
    scarry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        a0 = _mm512_load_epi32(s->data + i * VECLEN);
        b0 = _mm512_load_epi32(n->data + i * VECLEN);
        a0 = _mm512_mask_sbb_epi32(a0, scarry2, scarry, b0, &scarry);
        _mm512_store_epi32(c->data + i * VECLEN, a0);
    }
    _mm512_store_epi32(c->data + i * VECLEN, _mm512_set1_epi32(0));

    c->size = NWORDS;
    return;
}

void vecmodexp(bignum *d, bignum *b, bignum *e, bignum *m,
    bignum *s, bignum *one, monty *mdata)
{
    // d = b^e mod m
    // all b and e vector elements can be different.
    // all m elements can be different
    // proceed in a left to right fashion.
    int i, bit = NWORDS * DIGITBITS - 1;
    int j;
    // the bit string length.  needs to divide the word size (e.g., 32);
    int k = get_winsize();
    int mask;
    int bstr[VECLEN];
    uint8_t done[256];

    bignum **g = mdata->g;
    bignum *w = mdata->mtmp4;

    mask = 0;
    for (j = 0; j < k; j++)
    {
        mask = (mask << 1) | 1;
    }

    for (j = 16; j < (1 << k); j++)
    {
        done[j] = 0;
    }

    // precomputation.  smallest window is 4 bits, so these
    // computations can be pulled out of the loop
    vecCopy(b, g[1]);
    vecsqrmod_ptr(g[1], g[2], m, s, mdata);
    vecsqrmod_ptr(g[2], g[4], m, s, mdata);
    vecmulmod_ptr(g[2], b, g[3], m, s, mdata);
    vecsqrmod_ptr(g[3], g[6], m, s, mdata);
    vecmulmod_ptr(g[6], b, g[7], m, s, mdata);
    vecsqrmod_ptr(g[6], g[12], m, s, mdata);
    vecmulmod_ptr(g[4], b, g[5], m, s, mdata);
    vecsqrmod_ptr(g[4], g[8], m, s, mdata);
    vecsqrmod_ptr(g[5], g[10], m, s, mdata);
    vecsqrmod_ptr(g[7], g[14], m, s, mdata);
    vecmulmod_ptr(g[8], b, g[9], m, s, mdata);
    vecmulmod_ptr(g[10], b, g[11], m, s, mdata);
    vecmulmod_ptr(g[12], b, g[13], m, s, mdata);
    vecmulmod_ptr(g[14], b, g[15], m, s, mdata);
    //nsqr = 7;
    //nmul = 7;

    for (i = 16; i < (1 << k); i++)
    {
        int q = i;

        if (done[i])
            continue;

        vecmulmod_ptr(g[i - 1], b, g[i], m, s, mdata);

        while ((q * 2) < (1 << k))
        {
            if (!done[2 * q])
            {
                vecsqrmod_ptr(g[q], g[2 * q], m, s, mdata);
                done[2 * q] = 1;
            }
            q *= 2;
        }
    }

    vecCopyn(one, d, NWORDS);

    while (bit >= 0)
    {
        if (bit < k)
        {
            mask = 0x0;
            for (j = 0; j < (bit + 1); j++)
            {
                vecsqrmod_ptr(d, d, m, s, mdata);
                mask = (mask << 1) | 1;
            }

            for (j = 0; j < VECLEN; j++)
            {
                bstr[j] = (e->data[j + 0 * VECLEN]) & mask;
            }
        }
        else
        {
            for (j = 0; j < VECLEN; j++)
            {
                bstr[j] = get_bitwin(e, bit, k, j, mask);
            }

            if (bit < (DIGITBITS * NWORDS - 1))
            {
                for (j = 0; j < k; j++)
                {
                    vecsqrmod_ptr(d, d, m, s, mdata);
                }
            }
        }

        for (j = 0; j < VECLEN; j++)
        {
            copy_vec_lane(g[bstr[j]], w, j, NWORDS);
        }

        vecmulmod_ptr(d, w, d, m, s, mdata);

        bit -= k;
    }

    return;
}

uint32_t vec_mask_gte(uint32_t mask, bignum* u, bignum* v)
{
    // decide if each of the bignums in vec 'u' is >=
    // the corresponding bignum in vec 'v'.
    // return a mask of results.
    int i;
    __mmask16 mdecided = mask;
    __mmask16 mgte = 0;

    for (i = NWORDS - 1; i >= 0; --i)
    {
        __m512i a = _mm512_load_epi32(u->data + i * VECLEN);
        __m512i b = _mm512_load_epi32(v->data + i * VECLEN);

        mgte |= _mm512_mask_cmp_epu32_mask(~mdecided, a, b, _MM_CMPINT_GT);
        mdecided = mdecided | _mm512_mask_cmp_epu32_mask(~mdecided, a, b, _MM_CMPINT_LT);

        if (mdecided == 0xffff)
            break;
    }

    //equal if still undecided
    mgte |= ~mdecided;

    return (uint32_t)mgte;
}

uint32_t vec_gte(bignum * u, bignum * v)
{
    // decide if each of the bignums in vec 'u' is >=
    // the corresponding bignum in vec 'v'.
    // return a mask of results.
    int i;
    __mmask16 mdecided = 0;
    __mmask16 mgte = 0;

    for (i = NWORDS - 1; i >= 0; --i)
    {
        __m512i a = _mm512_load_epi32(u->data + i * VECLEN);
        __m512i b = _mm512_load_epi32(v->data + i * VECLEN);

        mgte |= _mm512_mask_cmp_epu32_mask(~mdecided, a, b, _MM_CMPINT_GT);
        mdecided = mdecided | _mm512_mask_cmp_epu32_mask(~mdecided, a, b, _MM_CMPINT_LT);

        if (mdecided == 0xffff)
            break;
    }

    //equal if still undecided
    mgte |= ~mdecided;

    return (uint32_t)mgte;
}

uint32_t vec_eq(base_t * u, base_t * v, int sz)
{
    // decide if each of the bignums in vec 'u' is >=
    // the corresponding bignum in vec 'v'.
    // return a mask of results.
    int i;
    __mmask16 meq = 0xffff;

    for (i = sz - 1; i >= 0; --i)
    {
        __m512i a = _mm512_load_epi32(u + i * VECLEN);
        __m512i b = _mm512_load_epi32(v + i * VECLEN);

        meq = _mm512_mask_cmp_epu32_mask(meq, a, b, _MM_CMPINT_EQ);

        if (meq == 0)
            break;
    }

    return (uint32_t)meq;
}

uint32_t vec_bignum_mask_lshift_1(bignum * u, uint32_t wmask)
{
    // return the left shift of bignum u by 1
    int i;
    __m512i nextcarry;
    __m512i carry = _mm512_setzero_epi32();
    __m512i word;

    for (i = 0; i < NWORDS; i++)
    {
        word = _mm512_load_epi32(u->data + i * VECLEN);
        // _mm512_and_epi32(word, highmask)  // not necessary to mask as the shift zero extends.
        nextcarry = _mm512_srli_epi32(word, 31);
        _mm512_mask_store_epi32(u->data + i * VECLEN, (__mmask16)wmask,
            _mm512_or_epi32(_mm512_slli_epi32(word, 1), carry));
        carry = nextcarry;
    }

    _mm512_mask_store_epi32(u->data + i * VECLEN, (__mmask16)wmask,
        _mm512_or_epi32(_mm512_slli_epi32(word, 1), carry));

    // return an overflow mask
    return wmask & _mm512_cmp_epi32_mask(carry, _mm512_setzero_epi32(), _MM_CMPINT_GT);
}

void vec_bignum_mask_rshift_1(bignum * u, uint32_t wmask)
{
    // return the right shift of bignum u by 1
    int i;
    __m512i nextcarry;
    __m512i carry = _mm512_setzero_epi32();
    __m512i lowmask = _mm512_set1_epi32(0x00000001);
    __m512i word;

    //carry = 0;
    //for (i = sb - 1; i >= 0; --i)
    //{
    //    nextcarry = (b->data[i] & mask) << y;
    //    a->data[i] = b->data[i] >> x | carry;
    //    carry = nextcarry;
    //}

    for (i = NWORDS - 1; i >= 0; i--)
    {
        word = _mm512_load_epi32(u->data + i * VECLEN);
        nextcarry = _mm512_slli_epi32(_mm512_and_epi32(word, lowmask), 31);
        _mm512_mask_store_epi32(u->data + i * VECLEN, (__mmask16)wmask,
            _mm512_or_epi32(_mm512_srli_epi32(word, 1), carry));
        carry = nextcarry;
    }

    return;
}

void vec_bignum_mask_sub(bignum *a, bignum *b, bignum *c, uint32_t wmask)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // s1 is of length VECLEN
    // a, b, c, and s1 are aligned
    // a and b are both positive
    // a >= b
    int i;

    __mmask16 carry = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;

    if (wmask == 0)
        return;

    // subtract the selected elements ('1' in the mask)
    carry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        avec = _mm512_load_epi32(a->data + i * VECLEN);
        bvec = _mm512_load_epi32(b->data + i * VECLEN);
        cvec = _mm512_sbb_epi32(avec, carry, bvec, &carry);
        _mm512_mask_store_epi32(c->data + i * VECLEN, (__mmask16)wmask, cvec);
    }

    if (carry)
    {
        // subtract any final borrows that exceed the size of b.
        _mm512_mask_store_epi32(c->data + i * VECLEN, (__mmask16)wmask & carry, _mm512_setzero_epi32());
    }

    return;
}

int vec_montgomery_setup(bignum * a, bignum *r, bignum *rhat, base_t *rho)
{
    __m512i x, b;
    __m512i four = _mm512_set1_epi32(4);
    __m512i two = _mm512_set1_epi32(2);
    int i, k, bits;
    uint32_t m1, m2;
    uint32_t barray[VECLEN];
    // fast inversion mod 2**k assuming odd integer input 'a'
    // 
    // Based on the fact that
    // 
    // XA = 1 (mod 2**n)  =>  (X(2-XA)) A = 1 (mod 2**2n)
    //                    =>  2*X*A - X*X*A*A = 1
    //                    =>  2*(1) - (1)     = 1
    // 
    b = _mm512_load_epi32(a->data);

    x = _mm512_add_epi32(b, two);
    x = _mm512_and_epi32(x, four);
    x = _mm512_slli_epi32(x, 1);
    x = _mm512_add_epi32(b, x);
    x = _mm512_mullo_epi32(x, _mm512_sub_epi32(two, _mm512_mullo_epi32(b, x)));
    x = _mm512_mullo_epi32(x, _mm512_sub_epi32(two, _mm512_mullo_epi32(b, x)));
    x = _mm512_mullo_epi32(x, _mm512_sub_epi32(two, _mm512_mullo_epi32(b, x)));
    _mm512_store_epi32(rho, _mm512_sub_epi32(_mm512_setzero_epi32(), x));

    // fixed precision... work with R = 2^(NWORDS*DIGITBITS)    
    vecClear(r);
    for (i = 0; i < VECLEN; i++)
        r->data[NWORDS * VECLEN + i] = 1;
    r->size = NWORDS + 1;

    // now compute r^2 % n so that we can do fast 
    // conversions into montgomery representation.
    bits = NWORDS * DIGITBITS - 1;

    // rhat initially starts out bigger than NWORDS, and the
    // rshift routine assumes numbers less than NWORDS, so
    // we do this first iteration outside the loop.
    vecClear(rhat);
    for (i = 0; i < VECLEN; i++)
        rhat->data[(NWORDS - 1) * VECLEN + i] = 0x80000000;
    bits += 2;
    rhat->size = NWORDS;

    for (i = 0; i < VECLEN; i++)
        barray[i] = bits;

    m1 = vec_gte(rhat, a);

    k = 0;
    while (m1 > 0)
    {
        for (i = 0; i < VECLEN; i++)
        {
            if (m1 & (1 << i))
                barray[i]++;
        }

        vec_bignum_mask_rshift_1(rhat, m1);
        m1 = vec_gte(rhat, a);
        bits++;
    }

    for (k = 0; k < bits; k++) {

        m1 = 0;
        for (i = 0; i < VECLEN; i++)
        {
            if (barray[i] > 0)
            {
                m1 |= (1 << i);
                barray[i]--;
            }
        }

        m2 = vec_bignum_mask_lshift_1(rhat, m1);
        m1 = vec_gte(rhat, a) | m2;
        vec_bignum_mask_sub(rhat, a, rhat, m1);
    }

    rhat->size = a->size;

    return 0;
}

