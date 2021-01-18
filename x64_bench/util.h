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

#ifndef _UTIL_H
#define _UTIL_H

// ============================================================================
// some standard headers
// ============================================================================
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>		// for uint32_t, etc.
#include <time.h>
#include <string.h>
#include <sys/types.h>
#if defined(WIN32)
#define WIN32_LEAN_AND_MEAN
#include <intrin.h>	
#include <malloc.h>
#include <windows.h>
#include <process.h>
#endif

#ifndef _MSC_VER
#include <sys/time.h>	//for gettimeofday using gcc
#include <unistd.h>
#endif


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
#define HBAR 6.58211928000e-7        // (ev * ns)
#define INV_HBAR 1.519267514702347e+06	// (eV * ns)^-1
#define INV_HBAR_RESIDUE 5.873436907300274 // 1/hbar mod 2*pi

#define INLINE __inline
#define LOWER(x) ((x) & HALFMASK)
#define UPPER(x) ((x) >> HALFBITS)
#define strto_uint64 strtoull
#define DEC 10
#define HEX 16
#define DEFINED 1
#ifdef NOTDEF
#undef NOTDEF
#endif

// portable 64-bit formatting and aligned memory 
#if defined(_MSC_VER) || defined(__MINGW32__)
#define PRId64 "I64d"
#define PRIu64 "I64u"
#define PRIx64 "I64x"

#define align_free _aligned_free
#define ALIGNED_MEM __declspec(align(64))     

#elif defined(__x86_64__)

#define align_free free
#if defined (__INTEL_COMPILER)
#define ALIGNED_MEM __declspec(align(64))
#else
#define ALIGNED_MEM __attribute__((aligned(64)))
#endif

#define PRId64 "ld"
#define PRIu64 "lu"
#define PRIx64 "lx"
#define BSCu "lu"
#define BSCx "lx"
#define BSCu0 "019lu"	// base string conversion with leading zeros
#define BSCx0 "019lx"	// base string conversion with leading zeros
#elif defined(__i386__)

#define align_free free
#if defined (__INTEL_COMPILER)
#define ALIGNED_MEM __declspec(align(64))
#else
#define ALIGNED_MEM __attribute__((aligned(64)))
#endif

#define PRId64 "lld"
#define PRIu64 "llu"
#define PRIx64 "llx"
#define BSCu "u"
#define BSCx "x"
#define BSCu0 "09u"
#define BSCx0 "09x"
#endif

#ifdef _MSC_VER
#define strto_uint64 _strtoui64
#else
#define strto_uint64 strtoull
#endif


// ============================================================================
// memory allocation
// ============================================================================
static __inline void * xmalloc_align(size_t len)
{
#if defined (_MSC_VER) || defined(__MINGW32__)
    void *ptr = _aligned_malloc(len, 64);
#elif defined (__APPLE__)
    void *ptr = malloc(len);
#elif defined (__GNUC__)
    void *ptr = memalign(64, len);
#define align_free free
#else
    void *ptr = malloc(len);
#endif

    if (ptr == NULL) {
        printf("failed to allocate %u aligned bytes\n", (uint32_t)len);
        exit(-1);
    }

    return ptr;
}

static __inline void * xmalloc(size_t len) {
    void *ptr = malloc(len);
    if (ptr == NULL) {
        printf("failed to allocate %u bytes\n", (uint32_t)len);
        exit(-1);
    }
    return ptr;
}

static __inline void * xcalloc(size_t num, size_t len) {
    void *ptr = calloc(num, len);
    if (ptr == NULL) {
        printf("failed to calloc %u bytes\n", (uint32_t)(num * len));
        exit(-1);
    }
    return ptr;
}

static __inline void * xrealloc(void *iptr, size_t len) {
    void *ptr = realloc(iptr, len);
    if (ptr == NULL) {
        printf("failed to reallocate %u bytes\n", (uint32_t)len);
        exit(-1);
    }
    return ptr;
}

// ============================================================================
// randomness
// ============================================================================
typedef struct
{
    uint32_t hi;
    uint32_t low;
} rand_t;

uint32_t spRand(uint64_t *state, uint32_t lower, uint32_t upper);
uint64_t spRand64(uint64_t *state);
uint64_t spRand64_range(uint64_t *state, uint64_t lower, uint64_t upper);
void get_random_seeds(rand_t *r);

rand_t g_rand;
uint64_t LCGSTATE;

// ============================================================================
// hashing
// ============================================================================
uint64_t hash64(uint64_t in);


// ============================================================================
// sorting (qsort)
// ============================================================================

static int qcomp_uint32(const void *x, const void *y)
{
    uint32_t *xx = (uint32_t *)x;
    uint32_t *yy = (uint32_t *)y;

    if (*xx > *yy)
        return 1;
    else if (*xx == *yy)
        return 0;
    else
        return -1;
}


// ============================================================================
// precision time
// ============================================================================
uint64_t read_clock(void);
uint64_t measure_processor_speed(int millisec);

#ifdef _MSC_VER
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64

struct timeval
{
    long tv_sec;
    long tv_usec;
};

struct timezone
{
    int  tz_minuteswest; /* minutes W of Greenwich */
    int  tz_dsttime;     /* type of dst correction */
};
#endif

double my_difftime(struct timeval *, struct timeval *);
#if defined (_MSC_VER)
int gettimeofday(struct timeval *tv, struct timezone *tz);

static void usleep(uint32_t usec)
{
    Sleep(usec / 1000);
}
#endif

#endif
