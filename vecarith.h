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
#define base_t uint32_t
#define DIGITBITS 32
#define VECLEN 16
#ifndef MAXBITS
#define MAXBITS 512
#endif
#define NWORDS (MAXBITS / DIGITBITS)
#define HALFBITS 16
#define HALFMASK 0xffff
#define MAXDIGIT 0xffffffffULL
#define HIBITMASK 0x80000000ULL
#define MAX_WINSIZE 8

// ============================================================================
// useful definitions
// ============================================================================
#define MIN(a,b) ((a) < (b)? (a) : (b))
#define MAX(a,b) ((a) > (b)? (a) : (b))
#define SIGN(a) ((a) < 0 ? -1 : 1)

#define INV_2_POW_48 3.5527136788005009293556213378906e-15
#define INV_2_POW_52 2.2204460492503130808472633361816e-16
#define INV_2_POW_64 5.4210108624275221700372640043497e-20
#define INV_2_POW_26 1.490116119384765625e-8
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
int vec_montgomery_setup(bignum * a, bignum *r, bignum *rhat, base_t *rho);
int get_winsize(void);
void montvecadd(bignum *a, bignum *b, bignum *c, bignum *n);
void montvecsub(bignum *a, bignum *b, bignum *c, bignum *n);
void vecmulmod(bignum *a, bignum *b, bignum *c, bignum *n, bignum *s, monty *mdata);
void vecsqrmod(bignum *a, bignum *c, bignum *n, bignum *s, monty *mdata);
void vecmodexp(bignum *d, bignum *b, bignum *e, bignum *m,
    bignum *s, bignum *one, monty *mdata);

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
uint32_t vec_eq(base_t * u, base_t * v, int sz);
uint32_t vec_bignum_mask_lshift_1(bignum * u, uint32_t wmask);
void vec_bignum_mask_rshift_1(bignum * u, uint32_t wmask);
void vec_bignum_mask_sub(bignum *a, bignum *b, bignum *c, uint32_t wmask);


