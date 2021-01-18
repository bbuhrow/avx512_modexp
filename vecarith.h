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

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <malloc.h>
#include <sys/time.h>	//for gettimeofday using gcc
#include <time.h>
#include <immintrin.h>

// ============================================================================
// vecarith config
// ============================================================================
#define DIGITBITS 32
#define base_t uint32_t
#define base_signed_t int32_t
#define HALFBITS 16
#define HALFMASK 0xffff
#define MAXDIGIT 0xffffffff
#define HIBITMASK 0x80000000
#define VECLEN 16

#ifndef MAXBITS
#define MAXBITS 512
#endif
#define NWORDS (MAXBITS / DIGITBITS)
#define DEFINED 1
#define MAX_WINSIZE 8

// ============================================================================
// useful definitions
// ============================================================================
#define MIN(a,b) ((a) < (b)? (a) : (b))
#define MAX(a,b) ((a) > (b)? (a) : (b))
#define SIGN(a) ((a) < 0 ? -1 : 1)

#define INV_2_POW_48 3.5527136788005009293556213378906e-15
#define INV_2_POW_64 5.4210108624275221700372640043497e-20
#define INV_2_POW_32 2.3283064365386962890625e-10
#define PI 3.1415926535897932384626433832795
#define LN2 0.69314718055994530941723212145818
#ifdef _MSC_VER
#define strto_uint64 _strtoui64
#else
#define strto_uint64 strtoull
#endif

// portable 64-bit formatting
#if defined(_MSC_VER) || defined(__MINGW32__)
#define PRId64 "I64d"
#define PRIu64 "I64u"
#define PRIx64 "I64x"
#elif defined(__x86_64__)
#define PRId64 "ld"
#define PRIu64 "lu"
#define PRIx64 "lx"
#define BSCu "lu"
#define BSCx "lx"
#define BSCu0 "019lu"	// base string conversion with leading zeros
#define BSCx0 "019lx"	// base string conversion with leading zeros
#elif defined(__i386__)
#define PRId64 "lld"
#define PRIu64 "llu"
#define PRIx64 "llx"
#define BSCu "u"
#define BSCx "x"
#define BSCu0 "09u"
#define BSCx0 "09x"
#endif


#if defined (__INTEL_COMPILER)
#define ALIGNED_MEM __declspec(align(64))
#else
#define ALIGNED_MEM __attribute__((aligned(64)))
#endif


// ============================================================================
// memory allocation
// ============================================================================
static __inline void * xmalloc_align(size_t len)
{
#if defined (_MSC_VER) || defined(__MINGW32__)
    void *ptr = _aligned_malloc(len, 64);
#define align_free _aligned_free
#elif defined (__APPLE__)
    void *ptr = malloc(len);
#elif defined (__GNUC__)
    void *ptr = memalign(64, len);
#define align_free free
#else
    void *ptr = malloc(len);
#endif

    if (ptr == NULL) {
        printf("failed to allocate %u aligned bytes\n", (uint32_t)len); fflush(stdout);
        exit(-1);
    }

    return ptr;
}

static __inline void * xmalloc(size_t len) {
    void *ptr = malloc(len);
    if (ptr == NULL) {
        printf("failed to allocate %u bytes\n", (uint32_t)len); fflush(stdout);
        exit(-1);
    }
    return ptr;
}

static __inline void * xcalloc(size_t num, size_t len) {
    void *ptr = calloc(num, len);
    if (ptr == NULL) {
        printf("failed to calloc %u bytes\n", (uint32_t)(num * len)); fflush(stdout);
        exit(-1);
    }
    return ptr;
}

static __inline void * xrealloc(void *iptr, size_t len) {
    void *ptr = realloc(iptr, len);
    if (ptr == NULL) {
        printf("failed to reallocate %u bytes\n", (uint32_t)len); fflush(stdout);
        exit(-1);
    }
    return ptr;
}

// ============================================================================
// vector bignum structure
// ============================================================================
typedef struct
{
    base_t *data;
    int size;
} bignum;

// ============================================================================
// montgomery arithmetic
// ============================================================================
typedef struct
{
    bignum *r;
    bignum *n;
    bignum *nhat;
    bignum *vnhat;
    bignum *rhat;
    bignum *rmask;
    bignum *one;
    bignum *mtmp1;
    bignum *mtmp2;
    bignum *mtmp3;
    bignum *mtmp4;
    bignum **g;             // storage for windowed method precomputation
    base_t *vrho;
    base_t rho;
} monty;

monty* monty_alloc(void);
void monty_free(monty *mdata); 
void monty_init_vec(monty *mdata, bignum * n, int verbose);
int get_winsize(void);
// 32-bit words, 16x
int vec_montgomery_setup(bignum * a, bignum *r, bignum *rhat, base_t *rho);
void vecmulmod(bignum *a, bignum *b, bignum *c, bignum *n, bignum *s, monty *mdata);
void vecsqrmod(bignum *a, bignum *c, bignum *n, bignum *s, monty *mdata);
void vecmodexp(bignum *d, bignum *b, bignum *e, bignum *m,
    bignum *s, bignum *one, monty *mdata);

void(*vecmulmod_ptr)(bignum *, bignum *, bignum *, bignum *, bignum *, monty *);
void(*vecsqrmod_ptr)(bignum *, bignum *, bignum *, bignum *, monty *);
int(*montsetup_ptr)(bignum *, bignum *, bignum *, base_t *);
void(*vecmodexp_ptr)(bignum *, bignum *, bignum *, bignum *, bignum *, bignum *, monty *m);

// ============================================================================
// vector bignum arithmetic and conversions
// ============================================================================
bignum * vecInit(void);
void vecCopy(bignum * src, bignum * dest);
void vecCopyn(bignum * src, bignum * dest, int size);
void vecClear(bignum *n);
void vecFree(bignum *n);
void broadcast_bignum_to_vec(bignum *src, bignum *vec_dest);
bignum * bignums_to_vec(bignum **src, int num);
void insert_bignum_in_vec(bignum *src, bignum *vec_dest, int num);
void extract_bignum_from_vec(bignum *vec_src, bignum *dest, int num);
void copy_vec_lane(bignum *src, bignum *dest, int num, int size);
uint32_t vec_gte(bignum * u, bignum * v);
uint32_t vec_mask_gte(uint32_t mask, bignum* u, bignum* v);
uint32_t vec_eq(base_t * u, base_t * v, int sz);
uint32_t vec_bignum_mask_lshift_1(bignum * u, uint32_t wmask);
void vec_bignum_mask_rshift_1(bignum * u, uint32_t wmask);
void vec_bignum_mask_sub(bignum *a, bignum *b, bignum *c, uint32_t wmask);

// ---------------------------------------------------------------------
// emulated instructions
// ---------------------------------------------------------------------
__m512i __inline _mm512_mulhi_epu32(__m512i a, __m512i b)
{
    __m512i t1 = _mm512_shuffle_epi32(a, 0xB1);
    __m512i t2 = _mm512_shuffle_epi32(b, 0xB1);
    __m512i evens = _mm512_mul_epu32(a, b);
    __m512i odds = _mm512_mul_epu32(t1, t2);
    //return _mm512_mask_mov_epi32(_mm512_shuffle_epi32(evens, 0xB1), 0xaaaa, odds);
    return _mm512_mask_mov_epi32(odds, 0x5555, _mm512_shuffle_epi32(evens, 0xB1));
}

__m512i __inline _mm512_mask_adc_epi32(__m512i a, __mmask16 m, __mmask16 c, __m512i b, __mmask16 *cout)
{
    __m512i t = _mm512_add_epi32(a, b);
    *cout = _mm512_cmplt_epu32_mask(t, a);
    __m512i t2 = _mm512_mask_add_epi32(a, m, t, _mm512_maskz_set1_epi32(c, 1));
    *cout = _mm512_kor(*cout, _mm512_mask_cmplt_epu32_mask(m, t2, t));
    return t2;
}

__m512i __inline _mm512_adc_epi32_test1(__m512i a, __mmask16 c, __m512i b, __mmask16 *cout)
{
    __m512i t = _mm512_add_epi32(a, b);
    *cout = _mm512_cmplt_epu32_mask(t, a);
    __m512i t2 = _mm512_add_epi32(t, _mm512_maskz_set1_epi32(c, 1));
    *cout = _mm512_kor(*cout, _mm512_cmplt_epu32_mask(t2, t));
    return t2;
}

__m512i __inline _mm512_adc_epi32_test2(__m512i a, __mmask16 c, __m512i b, __mmask16 *cout)
{
    // looks like a slightly improved data dependency chain... 
    // but it tested slower for 1024-b inputs...
    __m512i t = _mm512_add_epi32(a, b);
    __mmask16 gt0 = _mm512_kor(_mm512_test_epi32_mask(b, b), c);

    t = _mm512_add_epi32(t, _mm512_maskz_set1_epi32(c, 1));
    *cout = _mm512_kand(_mm512_cmple_epu32_mask(t, a), gt0);
    return t;
}

__m512i __inline _mm512_adc_epi32(__m512i a, __mmask16 c, __m512i b, __mmask16 *cout)
{
    __m512i t = _mm512_add_epi32(a, b);
    t = _mm512_add_epi32(t, _mm512_maskz_set1_epi32(c, 1));
    *cout = _mm512_cmplt_epu32_mask(t, a) | (_mm512_cmpeq_epu32_mask(t, a) & c);
    return t;
}

__m512i __inline _mm512_addcarry_epi32(__m512i a, __mmask16 c, __mmask16 *cout)
{
    __m512i t = _mm512_add_epi32(a, _mm512_maskz_set1_epi32(c, 1));
    *cout = _mm512_cmplt_epu32_mask(t, a);
    return t;
}

__m512i __inline _mm512_subborrow_epi32(__m512i a, __mmask16 c, __mmask16 *cout)
{
    __m512i t = _mm512_sub_epi32(a, _mm512_maskz_set1_epi32(c, 1));
    *cout = _mm512_cmpeq_epu32_mask(a, _mm512_setzero_epi32());
    return t;
}

__m512i __inline _mm512_mask_sbb_epi32(__m512i a, __mmask16 m, __mmask16 c, __m512i b, __mmask16 *cout)
{
    __m512i t = _mm512_sub_epi32(a, b);
    *cout = _mm512_mask_cmpgt_epu32_mask(m, t, a);
    __m512i t2 = _mm512_mask_sub_epi32(a, m, t, _mm512_maskz_set1_epi32(c, 1));
    *cout = _mm512_kor(*cout, _mm512_mask_cmpgt_epu32_mask(m, t2, t));
    return t2;
}

__m512i __inline _mm512_sbb_epi32(__m512i a, __mmask16 c, __m512i b, __mmask16 *cout)
{
    __m512i t = _mm512_sub_epi32(a, b);
    *cout = _mm512_cmpgt_epu32_mask(t, a);
    __m512i t2 = _mm512_sub_epi32(t, _mm512_maskz_set1_epi32(c, 1));
    *cout = _mm512_kor(*cout, _mm512_cmpgt_epu32_mask(t2, t));
    return t2;
}

__m512i __inline _mm512_sbb_epi64(__m512i a, __mmask8 c, __m512i b, __mmask8 *cout)
{
    __m512i t = _mm512_sub_epi64(a, b);
    *cout = _mm512_cmpgt_epu64_mask(t, a);
    __m512i t2 = _mm512_sub_epi64(t, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_kor(*cout, _mm512_cmpgt_epu64_mask(t2, t));
    return t2;
}

__m512i __inline _mm512_addsetc_epi32(__m512i a, __m512i b, __mmask16 *cout)
{
    __m512i t = _mm512_add_epi32(a, b);
    *cout = _mm512_cmplt_epu32_mask(t, a);
    return t;
}

__m512i __inline _mm512_subsetc_epi32(__m512i a, __m512i b, __mmask16 *cout)
{
    __m512i t = _mm512_sub_epi32(a, b);
    *cout = _mm512_cmpgt_epu32_mask(b, a);
    return t;
}

__inline void _mm512_epi32_to_eo64(__m512i a, __m512i *e64, __m512i *o64)
{
    *e64 = _mm512_maskz_mov_epi32(0x5555, a);
    *o64 = _mm512_maskz_mov_epi32(0x5555, _mm512_shuffle_epi32(a, 0xB1));
    return;
}

__inline __m512i _mm512_eo64lo_to_epi32(__m512i e64, __m512i o64)
{
    return _mm512_mask_blend_epi32(0xAAAA, e64, _mm512_shuffle_epi32(o64, 0xB1));
}

__inline __m512i _mm512_eo64hi_to_epi32(__m512i e64, __m512i o64)
{
    return _mm512_mask_blend_epi32(0xAAAA, _mm512_shuffle_epi32(e64, 0xB1), o64);
}

__inline void _mm512_mul_eo64_epi32(__m512i a, __m512i b, __m512i *e64, __m512i *o64)
{
    // multiply the 16-element 32-bit vectors a and b to produce two 8-element
    // 64-bit vector products e64 and o64, where e64 is the even elements
    // of a*b and o64 is the odd elements of a*b
    //__m512i t1 = _mm512_shuffle_epi32(a, 0xB1);
    //__m512i t2 = _mm512_shuffle_epi32(b, 0xB1);

    //_mm512_shuffle_epi32(a, 0xB1);
    //_mm512_shuffle_epi32(b, 0xB1);
    *e64 = _mm512_mul_epu32(a, b);
    *o64 = _mm512_mul_epu32(_mm512_shuffle_epi32(a, 0xB1), _mm512_shuffle_epi32(b, 0xB1));

    return;
}

#define _mm512_iseven_epi32(x) \
    _mm512_cmp_epi32_mask(_mm512_setzero_epi32(), _mm512_and_epi32((x), _mm512_set1_epi32(1)), _MM_CMPINT_EQ)

#define _mm512_isodd_epi32(x) \
    _mm512_cmp_epi32_mask(_mm512_set1_epi32(1), _mm512_and_epi32((x), _mm512_set1_epi32(1)), _MM_CMPINT_EQ)


