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

#include "monty_arith.h"
#include "x64_arith.h"

/* ==============================================================*/
/* Montgomery arithmetic setup */
/* ==============================================================*/

void to_monty(monty *mdata, bignum * x)
{
	// given a number x in normal (hexadecimal) representation, 
	// find its montgomery representation
       
	// this uses some precomputed monty constants
	// xhat = (x * r) mod n
    //  = x * R^2 / R mod n
    //  = REDC(x * R^2)
    mulmod_sos(mdata, x, mdata->rhat, mdata->mtmp2, mdata->mtmp1);
    if (zCompare(mdata->mtmp2, mdata->n))
    {
        zSub(mdata->mtmp2, mdata->n, x);
    }
    else
    {
        zCopy(mdata->mtmp2, x);
    }

	return;
}

monty * monty_alloc()
{
	//for a input modulus n, initialize constants for 
	//montogomery representation
	//this assumes that n is relatively prime to 2, i.e. is odd.	
	monty *mdata;
    int i;
	mdata = (monty *)malloc(sizeof(monty));

	// initialize space such that we can use vectorized
	mdata->n	 = zInit();
    mdata->np    = zInit();
	mdata->nhat  = zInit();
    mdata->vnhat = zInit();
	mdata->r     = zInit();
	mdata->one   = zInit();
	mdata->mtmp1 = zInit();
	mdata->mtmp2 = zInit();
	mdata->mtmp3 = zInit();
	mdata->mtmp4 = zInit();
    mdata->rmask = zInit();
	mdata->rhat = zInit();

	return mdata;
}

void monty_init(monty * mdata, bignum * n, int verbose)
{
    //for a input modulus n, initialize constants for 
    //montogomery representation
    //this assumes that n is relatively prime to 2, i.e. is odd.	
    uint64_t b = n->data[0];
    uint64_t x;
    int i, bits;

    if (verbose > 1)
        printf("initializing montgomery representation\n");

    zCopy(n, mdata->n);

    // invert (odd) n mod 2^64
    x = (((b + 2) & 4) << 1) + b; // here x*a==1 mod 2**4
    x *= 2 - b * x;               // here x*a==1 mod 2**8
    x *= 2 - b * x;               // here x*a==1 mod 2**16
    x *= 2 - b * x;               // here x*a==1 mod 2**32         
    x *= 2 - b * x;               // here x*a==1 mod 2**64

    mdata->rho = (uint64_t)((uint64_t)0 - ((uint64_t)x));

    // now compute r^2 % n so that we can do fast 
    // conversions into montgomery representation.
    bits = MAXBITS + 1;
    memset(mdata->rhat->data, 0, NWORDS * sizeof(uint64_t));
    mdata->rhat->data[NWORDS - 1] = 0x8000000000000000;
    mdata->rhat->size = NWORDS;
    while (zCompare(mdata->rhat, n) >= 0)
    {
        zShiftRight(mdata->rhat, mdata->rhat, 1);
        bits++;
    }

    for (i = 0; i < bits; i++) {
        zShiftLeft(mdata->rhat, mdata->rhat, 1);
        if (zCompare(mdata->rhat, n) >= 0)
            zSub(mdata->rhat, n, mdata->rhat);
    }
    mdata->rhat->data[NWORDS] = 0;

    zSet1(mdata->one, 1);
    to_monty(mdata, mdata->one);
    
    if (verbose > 1)
    {
        char *strn;

        strn = z2decstr(mdata->n);
        printf("N = %s\n", strn);
        free(strn);
        strn = z2decstr(mdata->rhat);
        printf("R^2 mod N = %s\n", strn);
        free(strn);
        strn = z2decstr(mdata->one);
        printf("ONE = %s\n", strn);
        free(strn);
        printf("rho = %016lx\n", mdata->rho);
    }

    return;
}

void monty_free(monty *mdata)
{
    int i;

	zFree(mdata->mtmp1);
	zFree(mdata->mtmp2);
	zFree(mdata->mtmp3);
	zFree(mdata->mtmp4);
	zFree(mdata->rhat);

	zFree(mdata->one);
	zFree(mdata->n);
    zFree(mdata->np);
	zFree(mdata->nhat);
	zFree(mdata->r);    

	return;
}

#define BLOCKSIZE 4
#define NBLOCKS (NWORDS / BLOCKSIZE)

// this memory will be lost
static uint64_t **col;
static int initialized = 0;

/* ==============================================================*/
/* BPS block subroutines */
/* ==============================================================*/

static void bps_fullblock_bsz4x2(uint64_t* A, uint64_t* B, uint64_t* C, uint64_t* M,
    int x, int y, uint64_t** s)
{
    // the bps multiplier can be improved by grouping a * b and s * n into
    // the same block accumulation (more reuse of the column accumulators).
    // Can't do this for squaring, because a * b is doubled but s * n isn't, and
    // the doubling happens for the entire column accumulator after the fullblocks
    // have been processed.
    int i;
    int aoff = (x - y) * BLOCKSIZE;
    int boff = (y - 1) * BLOCKSIZE + 1;

    for (i = 0; i < BLOCKSIZE; i++)
    {
        __asm__(
            "movq 0(%4), %%r8 \n\t"
            "movq 8(%4), %%r9 \n\t"
            "movq 16(%4), %%r10 \n\t"

            "movq %%r11, %%rax	\n\t"
            "mulq 0(%1) \n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq 16(%0), %%rax	\n\t"
            "mulq 8(%1)	\n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq 8(%0), %%rax	\n\t"
            "mulq 16(%1) \n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq 0(%0), %%rax	\n\t"
            "mulq 24(%1) \n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq 24(%2), %%rax	\n\t"
            "mulq 0(%3) \n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq 16(%2), %%rax	\n\t"
            "mulq 8(%3)	\n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq 8(%2), %%rax	\n\t"
            "mulq 16(%3) \n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq 0(%2), %%rax	\n\t"
            "mulq 24(%3) \n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq %%r8, 0(%4)		\n\t"
            "movq %%r9, 8(%4)		\n\t"
            "movq %%r10, 16(%4)		\n\t"
            :
        : "r"(&A[aoff]), "r"(&B[boff + i]), "r"(&C[aoff]), "r"(&M[boff + i]), "r"(s[i])
            : "rax", "rdx", "r8", "r9", "r10", "cc", "memory");
    }
    return;
}

static void bps_fullblock_bsz4(uint64_t* A, uint64_t* B, int x, int y, uint64_t** s)
{
    int i;
    int aoff = (x - y) * BLOCKSIZE;
    int boff = (y - 1) * BLOCKSIZE + 1;

    for (i = 0; i < BLOCKSIZE; i++)
    {
        __asm__(
            "movq 0(%2), %%r8 \n\t"
            "movq 8(%2), %%r9 \n\t"
            "movq 16(%2), %%r10 \n\t"

            "movq 24(%0), %%rax	\n\t"
            "mulq 0(%1) \n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq 16(%0), %%rax	\n\t"
            "mulq 8(%1)	\n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq 8(%0), %%rax	\n\t"
            "mulq 16(%1) \n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq 0(%0), %%rax	\n\t"
            "mulq 24(%1) \n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq %%r8, 0(%2)		\n\t"
            "movq %%r9, 8(%2)		\n\t"
            "movq %%r10, 16(%2)		\n\t"
            :
        : "r"(&A[aoff]), "r"(&B[boff + i]), "r"(s[i])
            : "rax", "rdx", "r8", "r9", "r10", "cc", "memory");
    }
    return;
}

static void bps_fullblock_bsz4x4(uint64_t* A, uint64_t* B, int x, int y, uint64_t** s)
{
    int i;
    int aoff = (x - y) * BLOCKSIZE;
    int boff = (y - 1) * BLOCKSIZE + 1;

    //for (i = 0; i < BLOCKSIZE; i++)
    {
        __asm__(
            /* Load all A inputs for this block, and the B inputs that are */
            /* used the most in the remaining available registers */
            "movq 24(%0), %%r11 \n\t"
            "movq 16(%0), %%r12 \n\t"
            "movq 8(%0), %%r13 \n\t"
            "movq 0(%0), %%r14 \n\t"
            "movq 24(%1), %%r15 \n\t"
            "movq 16(%1), %%rcx \n\t"

            /* i = 0 */
            "movq 0(%2), %%rbx \n\t"
            "movq 0(%%rbx), %%r8 \n\t"
            "movq 8(%%rbx), %%r9 \n\t"
            "movq 16(%%rbx), %%r10 \n\t"
            
            "movq %%r11, %%rax	\n\t"
            "mulq 0(%1) \n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq %%r12, %%rax	\n\t"
            "mulq 8(%1)	\n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq %%r13, %%rax	\n\t"
            "mulq %%rcx \n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq %%r14, %%rax	\n\t"
            "mulq %%r15 \n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq %%r8, 0(%%rbx)		\n\t"
            "movq %%r9, 8(%%rbx)		\n\t"
            "movq %%r10, 16(%%rbx)		\n\t"

            /* i = 1 */
            "movq 8(%2), %%rbx \n\t"
            "movq 0(%%rbx), %%r8 \n\t"
            "movq 8(%%rbx), %%r9 \n\t"
            "movq 16(%%rbx), %%r10 \n\t"

            "movq %%r11, %%rax	\n\t"
            "mulq 8(%1) \n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq %%r12, %%rax	\n\t"
            "mulq %%rcx	\n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq %%r13, %%rax	\n\t"
            "mulq %%r15 \n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq %%r14, %%rax	\n\t"
            "mulq 32(%1) \n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq %%r8, 0(%%rbx)		\n\t"
            "movq %%r9, 8(%%rbx)		\n\t"
            "movq %%r10, 16(%%rbx)		\n\t"

            /* i = 2 */
            "movq 16(%2), %%rbx \n\t"
            "movq 0(%%rbx), %%r8 \n\t"
            "movq 8(%%rbx), %%r9 \n\t"
            "movq 16(%%rbx), %%r10 \n\t"

            "movq %%r11, %%rax	\n\t"
            "mulq %%rcx \n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq %%r12, %%rax	\n\t"
            "mulq %%r15	\n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq %%r13, %%rax	\n\t"
            "mulq 32(%1) \n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq %%r14, %%rax	\n\t"
            "mulq 40(%1) \n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq %%r8, 0(%%rbx)		\n\t"
            "movq %%r9, 8(%%rbx)		\n\t"
            "movq %%r10, 16(%%rbx)		\n\t"

            /* i = 3 */
            "movq 24(%2), %%rbx \n\t"
            "movq 0(%%rbx), %%r8 \n\t"
            "movq 8(%%rbx), %%r9 \n\t"
            "movq 16(%%rbx), %%r10 \n\t"

            "movq %%r11, %%rax	\n\t"
            "mulq %%r15 \n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq %%r12, %%rax	\n\t"
            "mulq 32(%1)	\n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq %%r13, %%rax	\n\t"
            "mulq 40(%1) \n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq %%r14, %%rax	\n\t"
            "mulq 48(%1) \n\t"
            "addq %%rax, %%r8 \n\t"
            "adcq %%rdx, %%r9 \n\t"
            "adcq $0, %%r10 \n\t"

            "movq %%r8, 0(%%rbx)		\n\t"
            "movq %%r9, 8(%%rbx)		\n\t"
            "movq %%r10, 16(%%rbx)		\n\t"
            :
        : "r"(&A[aoff]), "r"(&B[boff]), "r"(s)
            : "rax", "rbx", "rcx", "rdx", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "cc", "memory");
    }
    return;
}

static void bps_fullblock_bsz4x2_b(uint64_t* A, uint64_t* B, int x, int y, uint64_t** s)
{
    int i;
    int aoff = (x - y) * BLOCKSIZE;
    int boff = (y - 1) * BLOCKSIZE + 1;

    // unrolling by 2 by itself does nothing (actually hurts).  But
    // if we also use mulx so we can use independent adcx adox,
    // maybe that will help.  
    // Results: no, this is even slower.
    for (i = 0; i < BLOCKSIZE; i+=2)
    {
        __asm__(
            "movq 0(%2), %%rbx \n\t"
            "movq 0(%%rbx), %%r8 \n\t"
            "movq 8(%%rbx), %%r9 \n\t"
            "movq 16(%%rbx), %%r10 \n\t"

            "movq 8(%2), %%rbx \n\t"
            "movq 0(%%rbx), %%r11 \n\t"
            "movq 8(%%rbx), %%r12 \n\t"
            "movq 16(%%rbx), %%r13 \n\t"

            "movq 24(%0), %%rdx	\n\t"
            "mulx 0(%1), %%r14, %%r15 \n\t"
            "mulx 8(%1), %%rbx, %%rcx \n\t"
            "xorq %%rdx, %%rdx \n\t"
            "adcx %%r14, %%r8 \n\t"
            "adox %%rbx, %%r11 \n\t"
            "adcx %%r15, %%r9 \n\t"
            "adox %%rcx, %%r12 \n\t"
            "adcx %%rdx, %%r10 \n\t"
            "adox %%rdx, %%r13 \n\t"

            "movq 16(%0), %%rdx	\n\t"
            "mulx 8(%1), %%r14, %%r15 \n\t"
            "mulx 16(%1), %%rbx, %%rcx \n\t"
            "xorq %%rdx, %%rdx \n\t"
            "adcx %%r14, %%r8 \n\t"
            "adox %%rbx, %%r11 \n\t"
            "adcx %%r15, %%r9 \n\t"
            "adox %%rcx, %%r12 \n\t"
            "adcx %%rdx, %%r10 \n\t"
            "adox %%rdx, %%r13 \n\t"

            "movq 8(%0), %%rdx	\n\t"
            "mulx 16(%1), %%r14, %%r15 \n\t"
            "mulx 24(%1), %%rbx, %%rcx \n\t"
            "xorq %%rdx, %%rdx \n\t"
            "adcx %%r14, %%r8 \n\t"
            "adox %%rbx, %%r11 \n\t"
            "adcx %%r15, %%r9 \n\t"
            "adox %%rcx, %%r12 \n\t"
            "adcx %%rdx, %%r10 \n\t"
            "adox %%rdx, %%r13 \n\t"

            "movq 0(%0), %%rdx	\n\t"
            "mulx 24(%1), %%r14, %%r15 \n\t"
            "mulx 32(%1), %%rbx, %%rcx \n\t"
            "xorq %%rdx, %%rdx \n\t"
            "adcx %%r14, %%r8 \n\t"
            "adox %%rbx, %%r11 \n\t"
            "adcx %%r15, %%r9 \n\t"
            "adox %%rcx, %%r12 \n\t"
            "adcx %%rdx, %%r10 \n\t"
            "adox %%rdx, %%r13 \n\t"

            "movq 0(%2), %%rbx \n\t"
            "movq %%r8, 0(%%rbx)		\n\t"
            "movq %%r9, 8(%%rbx)		\n\t"
            "movq %%r10, 16(%%rbx)		\n\t"

            "movq 8(%2), %%rbx \n\t"
            "movq %%r11, 0(%%rbx)		\n\t"
            "movq %%r12, 8(%%rbx)		\n\t"
            "movq %%r13, 16(%%rbx)		\n\t"
            :
        : "r"(&A[aoff]), "r"(&B[boff + i]), "r"(&s[i])
            : "rax", "rbx", "rcx", "rdx", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "cc", "memory");
    }
    return;
}

static void bps_pblock1p_bsz4(uint64_t* A, uint64_t* B, int x, uint64_t** s)
{
    // algorithm PartialBlock1', used on odd iterations, for a blocksize of 4.
    // The following terms will be doubled later, at the end of the column accumulation.
    // The square terms are handled separately, after the doubling step.

    //spMulAddc(A[x / 2 * BLOCKSIZE + 1], B[x / 2 * BLOCKSIZE + 3], s[0]);
    //spMulAddc(A[x / 2 * BLOCKSIZE + 0], B[x / 2 * BLOCKSIZE + 4], s[0]);

    __asm__(
        "movq 0(%2), %%r8 \n\t"
        "movq 8(%2), %%r9 \n\t"
        "movq 16(%2), %%r10 \n\t"

        "movq 8(%0), %%rax	\n\t"
        "mulq 24(%1) \n\t"
        "addq %%rax, %%r8 \n\t"
        "adcq %%rdx, %%r9 \n\t"
        "adcq $0, %%r10 \n\t"

        "movq 0(%0), %%rax	\n\t"
        "mulq 32(%1)	\n\t"
        "addq %%rax, %%r8 \n\t"
        "adcq %%rdx, %%r9 \n\t"
        "adcq $0, %%r10 \n\t"

        "movq %%r8, 0(%2)		\n\t"
        "movq %%r9, 8(%2)		\n\t"
        "movq %%r10, 16(%2)		\n\t"
        :
    : "r"(&A[x / 2 * BLOCKSIZE]), "r"(&B[x / 2 * BLOCKSIZE]), "r"(s[0])
        : "rax", "rdx", "r8", "r9", "r10", "cc", "memory");

    //spMulAddc(A[x / 2 * BLOCKSIZE + 2], B[x / 2 * BLOCKSIZE + 3], s[1]);
    //spMulAddc(A[x / 2 * BLOCKSIZE + 1], B[x / 2 * BLOCKSIZE + 4], s[1]);
    //spMulAddc(A[x / 2 * BLOCKSIZE + 0], B[x / 2 * BLOCKSIZE + 5], s[1]);

    __asm__(
        "movq 0(%2), %%r8 \n\t"
        "movq 8(%2), %%r9 \n\t"
        "movq 16(%2), %%r10 \n\t"

        "movq 16(%0), %%rax	\n\t"
        "mulq 24(%1) \n\t"
        "addq %%rax, %%r8 \n\t"
        "adcq %%rdx, %%r9 \n\t"
        "adcq $0, %%r10 \n\t"

        "movq 8(%0), %%rax	\n\t"
        "mulq 32(%1)	\n\t"
        "addq %%rax, %%r8 \n\t"
        "adcq %%rdx, %%r9 \n\t"
        "adcq $0, %%r10 \n\t"

        "movq 0(%0), %%rax	\n\t"
        "mulq 40(%1) \n\t"
        "addq %%rax, %%r8 \n\t"
        "adcq %%rdx, %%r9 \n\t"
        "adcq $0, %%r10 \n\t"

        "movq %%r8, 0(%2)		\n\t"
        "movq %%r9, 8(%2)		\n\t"
        "movq %%r10, 16(%2)		\n\t"
        :
    : "r"(&A[x / 2 * BLOCKSIZE]), "r"(&B[x / 2 * BLOCKSIZE]), "r"(s[1])
        : "rax", "rdx", "r8", "r9", "r10", "cc", "memory");

    //spMulAddc(A[x / 2 * BLOCKSIZE + 2], B[x / 2 * BLOCKSIZE + 4], s[2]);
    //spMulAddc(A[x / 2 * BLOCKSIZE + 1], B[x / 2 * BLOCKSIZE + 5], s[2]);
    //spMulAddc(A[x / 2 * BLOCKSIZE + 0], B[x / 2 * BLOCKSIZE + 6], s[2]);

    __asm__(
        "movq 0(%2), %%r8 \n\t"
        "movq 8(%2), %%r9 \n\t"
        "movq 16(%2), %%r10 \n\t"

        "movq 16(%0), %%rax	\n\t"
        "mulq 32(%1) \n\t"
        "addq %%rax, %%r8 \n\t"
        "adcq %%rdx, %%r9 \n\t"
        "adcq $0, %%r10 \n\t"

        "movq 8(%0), %%rax	\n\t"
        "mulq 40(%1)	\n\t"
        "addq %%rax, %%r8 \n\t"
        "adcq %%rdx, %%r9 \n\t"
        "adcq $0, %%r10 \n\t"

        "movq 0(%0), %%rax	\n\t"
        "mulq 48(%1) \n\t"
        "addq %%rax, %%r8 \n\t"
        "adcq %%rdx, %%r9 \n\t"
        "adcq $0, %%r10 \n\t"

        "movq %%r8, 0(%2)		\n\t"
        "movq %%r9, 8(%2)		\n\t"
        "movq %%r10, 16(%2)		\n\t"
        :
    : "r"(&A[x / 2 * BLOCKSIZE]), "r"(&B[x / 2 * BLOCKSIZE]), "r"(s[2])
        : "rax", "rdx", "r8", "r9", "r10", "cc", "memory");

    //spMulAddc(A[x / 2 * BLOCKSIZE + 1], B[x / 2 * BLOCKSIZE + 6], s[3]);
    //spMulAddc(A[x / 2 * BLOCKSIZE + 0], B[x / 2 * BLOCKSIZE + 7], s[3]);
    //spMulAddc(A[x / 2 * BLOCKSIZE + 3], B[x / 2 * BLOCKSIZE + 4], s[3]);
    //spMulAddc(A[x / 2 * BLOCKSIZE + 2], B[x / 2 * BLOCKSIZE + 5], s[3]);

    __asm__(
        "movq 0(%2), %%r8 \n\t"
        "movq 8(%2), %%r9 \n\t"
        "movq 16(%2), %%r10 \n\t"

        "movq 8(%0), %%rax	\n\t"
        "mulq 48(%1) \n\t"
        "addq %%rax, %%r8 \n\t"
        "adcq %%rdx, %%r9 \n\t"
        "adcq $0, %%r10 \n\t"

        "movq 0(%0), %%rax	\n\t"
        "mulq 56(%1)	\n\t"
        "addq %%rax, %%r8 \n\t"
        "adcq %%rdx, %%r9 \n\t"
        "adcq $0, %%r10 \n\t"

        "movq 24(%0), %%rax	\n\t"
        "mulq 32(%1) \n\t"
        "addq %%rax, %%r8 \n\t"
        "adcq %%rdx, %%r9 \n\t"
        "adcq $0, %%r10 \n\t"

        "movq 16(%0), %%rax	\n\t"
        "mulq 40(%1) \n\t"
        "addq %%rax, %%r8 \n\t"
        "adcq %%rdx, %%r9 \n\t"
        "adcq $0, %%r10 \n\t"

        "movq %%r8, 0(%2)		\n\t"
        "movq %%r9, 8(%2)		\n\t"
        "movq %%r10, 16(%2)		\n\t"
        :
    : "r"(&A[x / 2 * BLOCKSIZE]), "r"(&B[x / 2 * BLOCKSIZE]), "r"(s[3])
        : "rax", "rdx", "r8", "r9", "r10", "cc", "memory");

    return;
}

static void bps_pblock1pp_bsz4(uint64_t* A, uint64_t* B, int x, uint64_t** s)
{
    // algorithm PartialBlock1'', used on even iterations, for a blocksize of 4.
    // The following terms will be doubled later, at the end of the column accumulation.
    // The square terms are handled separately, after the doubling step.
    // Given the small number of terms, don't bother combining the block into
    // a single asm statement.
    spMulAddc(A[x / 2 * BLOCKSIZE + 0], B[x / 2 * BLOCKSIZE + 1], s[1]);
    spMulAddc(A[x / 2 * BLOCKSIZE + 0], B[x / 2 * BLOCKSIZE + 2], s[2]);
    spMulAddc(A[x / 2 * BLOCKSIZE + 1], B[x / 2 * BLOCKSIZE + 2], s[3]);
    spMulAddc(A[x / 2 * BLOCKSIZE + 0], B[x / 2 * BLOCKSIZE + 3], s[3]);

    return;
}

static void bps_pblock2pp_bsz4(uint64_t* A, uint64_t* B, int x, uint64_t** s)
{
    // algorithm PartialBlock2', used on even iterations, for a blocksize of 4.
    // The following terms will be doubled later, at the end of the column accumulation.
    // The square terms are handled separately, after the doubling step.

    // terms with pairs outside the block
    //spMulAddc(A[x / 2 * BLOCKSIZE + 1], B[x / 2 * BLOCKSIZE - 1], s[0]);
    //spMulAddc(A[x / 2 * BLOCKSIZE + 2], B[x / 2 * BLOCKSIZE - 2], s[0]);
    //spMulAddc(A[x / 2 * BLOCKSIZE + 3], B[x / 2 * BLOCKSIZE - 3], s[0]);

    __asm__(
        "movq 0(%2), %%r8 \n\t"
        "movq 8(%2), %%r9 \n\t"
        "movq 16(%2), %%r10 \n\t"

        "movq 8(%0), %%rax	\n\t"
        "mulq -8(%1) \n\t"
        "addq %%rax, %%r8 \n\t"
        "adcq %%rdx, %%r9 \n\t"
        "adcq $0, %%r10 \n\t"

        "movq 16(%0), %%rax	\n\t"
        "mulq -16(%1)	\n\t"
        "addq %%rax, %%r8 \n\t"
        "adcq %%rdx, %%r9 \n\t"
        "adcq $0, %%r10 \n\t"

        "movq 24(%0), %%rax	\n\t"
        "mulq -24(%1)	\n\t"
        "addq %%rax, %%r8 \n\t"
        "adcq %%rdx, %%r9 \n\t"
        "adcq $0, %%r10 \n\t"

        "movq %%r8, 0(%2)		\n\t"
        "movq %%r9, 8(%2)		\n\t"
        "movq %%r10, 16(%2)		\n\t"
        :
    : "r"(&A[x / 2 * BLOCKSIZE]), "r"(&B[x / 2 * BLOCKSIZE]), "r"(s[0])
        : "rax", "rdx", "r8", "r9", "r10", "cc", "memory");

    //spMulAddc(A[x / 2 * BLOCKSIZE + 2], B[x / 2 * BLOCKSIZE - 1], s[1]);
    //spMulAddc(A[x / 2 * BLOCKSIZE + 3], B[x / 2 * BLOCKSIZE - 2], s[1]);
    //spMulAddc(A[x / 2 * BLOCKSIZE + 1], B[x / 2 * BLOCKSIZE - 0], s[1]);

    __asm__(
        "movq 0(%2), %%r8 \n\t"
        "movq 8(%2), %%r9 \n\t"
        "movq 16(%2), %%r10 \n\t"

        "movq 16(%0), %%rax	\n\t"
        "mulq -8(%1) \n\t"
        "addq %%rax, %%r8 \n\t"
        "adcq %%rdx, %%r9 \n\t"
        "adcq $0, %%r10 \n\t"

        "movq 24(%0), %%rax	\n\t"
        "mulq -16(%1)	\n\t"
        "addq %%rax, %%r8 \n\t"
        "adcq %%rdx, %%r9 \n\t"
        "adcq $0, %%r10 \n\t"

        "movq 8(%0), %%rax	\n\t"
        "mulq 0(%1)	\n\t"
        "addq %%rax, %%r8 \n\t"
        "adcq %%rdx, %%r9 \n\t"
        "adcq $0, %%r10 \n\t"

        "movq %%r8, 0(%2)		\n\t"
        "movq %%r9, 8(%2)		\n\t"
        "movq %%r10, 16(%2)		\n\t"
        :
    : "r"(&A[x / 2 * BLOCKSIZE]), "r"(&B[x / 2 * BLOCKSIZE]), "r"(s[1])
        : "rax", "rdx", "r8", "r9", "r10", "cc", "memory");

    //spMulAddc(A[x / 2 * BLOCKSIZE + 2], B[x / 2 * BLOCKSIZE - 0], s[2]);
    //spMulAddc(A[x / 2 * BLOCKSIZE + 3], B[x / 2 * BLOCKSIZE - 1], s[2]);

    __asm__(
        "movq 0(%2), %%r8 \n\t"
        "movq 8(%2), %%r9 \n\t"
        "movq 16(%2), %%r10 \n\t"

        "movq 16(%0), %%rax	\n\t"
        "mulq 0(%1) \n\t"
        "addq %%rax, %%r8 \n\t"
        "adcq %%rdx, %%r9 \n\t"
        "adcq $0, %%r10 \n\t"

        "movq 24(%0), %%rax	\n\t"
        "mulq -8(%1)	\n\t"
        "addq %%rax, %%r8 \n\t"
        "adcq %%rdx, %%r9 \n\t"
        "adcq $0, %%r10 \n\t"

        "movq %%r8, 0(%2)		\n\t"
        "movq %%r9, 8(%2)		\n\t"
        "movq %%r10, 16(%2)		\n\t"
        :
    : "r"(&A[x / 2 * BLOCKSIZE]), "r"(&B[x / 2 * BLOCKSIZE]), "r"(s[2])
        : "rax", "rdx", "r8", "r9", "r10", "cc", "memory");

    //spMulAddc(A[x / 2 * BLOCKSIZE + 3], B[x / 2 * BLOCKSIZE - 0], s[3]);
    //spMulAddc(A[x / 2 * BLOCKSIZE + 2], B[x / 2 * BLOCKSIZE + 1], s[3]);

    __asm__(
        "movq 0(%2), %%r8 \n\t"
        "movq 8(%2), %%r9 \n\t"
        "movq 16(%2), %%r10 \n\t"

        "movq 24(%0), %%rax	\n\t"
        "mulq 0(%1) \n\t"
        "addq %%rax, %%r8 \n\t"
        "adcq %%rdx, %%r9 \n\t"
        "adcq $0, %%r10 \n\t"

        "movq 16(%0), %%rax	\n\t"
        "mulq 8(%1)	\n\t"
        "addq %%rax, %%r8 \n\t"
        "adcq %%rdx, %%r9 \n\t"
        "adcq $0, %%r10 \n\t"

        "movq %%r8, 0(%2)		\n\t"
        "movq %%r9, 8(%2)		\n\t"
        "movq %%r10, 16(%2)		\n\t"
        :
    : "r"(&A[x / 2 * BLOCKSIZE]), "r"(&B[x / 2 * BLOCKSIZE]), "r"(s[3])
        : "rax", "rdx", "r8", "r9", "r10", "cc", "memory");

    return;
}

static void bps_pblock2p_bsz4(uint64_t* A, uint64_t* B, int x, uint64_t** s)
{
    // algorithm PartialBlock2'', used on odd iterations, for a blocksize of 4.
    // The following terms will be doubled later, at the end of the column accumulation.
    // The square terms are handled separately, after the doubling step.
    // Given the small number of terms, don't bother combining the block into
    // a single asm statement.
    spMulAddc(A[x / 2 * BLOCKSIZE + 3], B[x / 2 * BLOCKSIZE + 1], s[0]);
    spMulAddc(A[x / 2 * BLOCKSIZE + 3], B[x / 2 * BLOCKSIZE + 2], s[1]);

    return;
}

static void bps_pblock1(uint64_t* A, uint64_t* B, int x, uint64_t** s)
{
    int i, j;
    uint64_t a[BLOCKSIZE];
    uint64_t b[BLOCKSIZE];

    for (i = 0; i < BLOCKSIZE; i++)
    {
        a[i] = A[x * BLOCKSIZE + i];
    }

    for (i = 0; i < BLOCKSIZE; i++)
    {
        b[i] = B[i];
    }

    // B=4, 10 it = n*(n+1)/2
    // B=5, 15 it
    // B=6, 21 it
    for (i = 0; i < BLOCKSIZE; i++)
    {
        for (j = 0; j <= i; j++)
        {
            // 5 reads, 3 writes, 3 adds, 1 mul
            spMulAddc(a[i - j], b[j], s[i]);
        }
    }

    return;
}

static void bps_pblock2(uint64_t* A, uint64_t* B, int x, uint64_t** s)
{
    int i, j;
    uint64_t a[BLOCKSIZE - 1];
    uint64_t b[BLOCKSIZE - 1];

    for (i = 1; i < BLOCKSIZE; i++)
    {
        a[i - 1] = A[(x - NBLOCKS) * BLOCKSIZE + i];
    }

    for (i = 1; i < BLOCKSIZE; i++)
    {
        b[i - 1] = B[NWORDS - i];
    }

    // B=4: 3+2+1 = 6 = B(B-1) iterations
    for (i = 0; i <= BLOCKSIZE - 2; i++)
    {
        for (j = i + 1; j < BLOCKSIZE; j++)
        {
            //s[i] = s[i] + a[j] * b[j - i];
            // 5 reads, 3 writes, 3 adds, 1 mul
            spMulAddc(a[j - 1], b[j - i - 1], s[i]);
        }
    }
    return;
}

static void bps_final(uint64_t* A, uint64_t* B, uint64_t** s, int i, uint64_t* acc, uint64_t rho)
{
    int p, q;

    // B*(B+1)/2 * (5 reads, 3 writes, 3 adds, 1 mul)
    // plus B * (7 reads, 3 adds, 1 mul, 4 writes)
    for (p = 0; p < BLOCKSIZE; p++)
    {
        //acc[0] += s[p];
        int c = 0;
        c = _addcarry_u64(c, acc[0], s[p][0], &acc[0]);
        c = _addcarry_u64(c, acc[1], s[p][1], &acc[1]);
        _addcarry_u64(c, acc[2], s[p][2], &acc[2]);

        for (q = 0; q < p; q++)
        {
            // 5 reads, 3 writes, 3 adds, 1 mul
            spMulAddc(A[i * BLOCKSIZE + q], B[p - q], acc);
        }

        // lower word only
        // 1 read, 1 mul, 1 write
        A[i * BLOCKSIZE + p] = acc[0] * rho;

        // multiply-accumulate with a final rotation of the accumulator
        // 5 reads, 3 writes, 3 adds, 1 mul
        spMulAddcr(A[i * BLOCKSIZE + p], B[0], acc);
    }

    return;
}


/* ==============================================================*/
/* Mulmod routines */
/* ==============================================================*/


void mulmod_sos(monty *mdata, bignum * u, bignum * v, bignum * w, bignum * s)
{
    // separate operand scanning. 
    // first multiply, then reduce
    int i, j;
    uint64_t c, p;

    memset(s->data, 0, (2 * NWORDS + 1) * sizeof(uint64_t));

    for (i = 0; i < NWORDS; i++)
    {
        c = 0;
        for (j = 0; j < NWORDS; j++)
        {
            spMulAdd2(u->data[j], v->data[i], s->data[i + j], c, &s->data[i + j], &c);
        }
        mpAdd1(s->data + i + j, c, s->data + i + j, 2 * NWORDS - i - j);
    }
    s->size = 2 * NWORDS;

    for (i = 0; i < NWORDS; i++)
    {
        spMul(s->data[i], mdata->rho, &p, &c);
        c = 0;
        for (j = 0; j < NWORDS; j++)
        {
            spMulAdd2(mdata->n->data[j], p, s->data[i + j], c, &s->data[i + j], &c);
        }
        mpAdd1(s->data + i + j, c, s->data + i + j, 2 * NWORDS - i - j);
    }
    s->size = 2 * NWORDS;

    // only get the result < 2^r, not strictly less than N
    if (s->data[2 * NWORDS])
    {
        mpSub(s->data + NWORDS, mdata->n->data, w->data, NWORDS);
    }
    else
    {
        memcpy(w->data, s->data + NWORDS, NWORDS * sizeof(uint64_t));
    }
    w->size = NWORDS;

    return;
}

void mulmod_cios(monty *mdata, bignum * u, bignum * v, bignum * w, bignum * s)
{
    // integrate multiply and reduction steps, alternating
    // between iterations of the outer loops.
    int i, j;
    uint64_t c, p;

    memset(s->data, 0, (NWORDS + 2) * sizeof(uint64_t));

    for (i = 0; i < NWORDS; i++)
    {
        c = 0;
        for (j = 0; j < NWORDS; j++)
        {
            spMulAdd2x(u->data[j], v->data[i], s->data[j], c, &s->data[j], &c);
        }
        s->data[j] += c;
        s->data[j + 1] = 0;

        p = s->data[0] * mdata->rho;

        spMulAdd1(mdata->n->data[0], p, s->data[0], &s->data[0], &c);
        for (j = 1; j < NWORDS; j++)
        {
            spMulAdd2x(mdata->n->data[j], p, s->data[j], c, &s->data[j - 1], &c);
        }
        spAdd(s->data[j], c, &s->data[j - 1], &s->data[j]);
        s->data[j + 1] = 0;
    }
    s->size = NWORDS;

    // only get the result < 2^r, not strictly less than N
    if (s->data[NWORDS])
    {
        mpSub(s->data, mdata->n->data, w->data, NWORDS);        
    }
    else
    {
        memcpy(w->data, s->data, NWORDS * sizeof(uint64_t));
    }
    w->data[NWORDS] = 0;
    w->data[NWORDS+1] = 0;
    w->size = NWORDS;

    return;
}

void mulmod_fios(monty *mdata, bignum * u, bignum * v, bignum * w, bignum * s)
{
    // integrate multiply and reduction steps, alternating
    // between iterations of the inner loops.
    int i, j;
    uint64_t c, c2, p, p2;

    memset(s->data, 0, (NWORDS + 2) * sizeof(uint64_t));

    for (i = 0; i < NWORDS; i++)
    {
        spMulAdd1(u->data[0], v->data[i], s->data[0], &s->data[0], &c);
        spAdd(s->data[1], c, s->data + 1, &c2);     // c2 applies to j+2 (next j+1)
        spMul(s->data[0], mdata->rho, &p, &c);
        spMulAdd1(mdata->n->data[0], p, s->data[0], &s->data[0], &c);

        for (j = 1; j < NWORDS; j++)
        {
            spMulAdd2(u->data[j], v->data[i], s->data[j], c, &p2, &c);
            spAdd3(s->data[j + 1], c, c2, s->data + j + 1, &c2);     // c2 applies to j+2 (next j+1)
            spMulAdd1(mdata->n->data[j], p, p2, &s->data[j - 1], &c);  // c applies to j+1 (next j)
        }
        spAdd(s->data[j], c, &s->data[j - 1], &s->data[j]);
        s->data[j + 1] = 0;
    }

    // only get the result < 2^r, not strictly less than N
    if (s->data[NWORDS])
    {
        mpSub(s->data, mdata->n->data, w->data, NWORDS);
    }
    else
    {
        memcpy(w->data, s->data, NWORDS * sizeof(uint64_t));
    }
    w->size = NWORDS;

    return;
}

void mulmod_fips(monty *mdata, bignum * u, bignum * v, bignum * w, bignum * s)
{
    // integrate multiply and reduction steps, alternating
    // between iterations of the inner loops while keeping
    // computations in product scanning form.
    int i, j;
    uint64_t c, p;
    uint64_t t[3] = { 0, 0, 0 };

    memset(s->data, 0, (NWORDS + 1) * sizeof(uint64_t));

    // muls
    // nwords = 4: 2*2*6+4*3 muls  = 36 = 2*s^2 + s = 36
    // nwords = 5: 2*2*10+5*3 muls = 55 = 2*s^2 + s = 55
    // nwords = 6: 2*2*15+6*3 muls = 78 = 2*s^2 + s = 78
    // fewer reads/writes than in Koc because we r/w the accumulator 
    // only once for a*b and s*n for a given i,j
    // nwords = 4: 7*2*6+4*15  reads = 144 = 7*s^2 + 8s
    // nwords = 5: 7*2*10+5*15 reads = 215 = 7*s^2 + 8s
    // nwords = 6: 7*2*15+6*15 reads = 300 = 7*s^2 + 8s
    // writes
    // nwords = 4: 3*2*6+4*11  writes = 80 = 3*s^2 + 8s
    // nwords = 5: 3*2*10+5*11 writes = 115 = 3*s^2 + 8s
    // nwords = 6: 3*2*15+6*11 writes = 156 = 3*s^2 + 8s
    // adds - we find fewer adds are required than 6*s^2 + 2s + 2, 
    // only 6*s^2.  maybe because of adc?
    // nwords = 4: 6*2*6+4*6  adds = 96  = 6*s^2 + 2s + 2 = 106
    // nwords = 5: 6*2*10+5*6 adds = 150 = 6*s^2 + 2s + 2 = 162
    // nwords = 6: 6*2*15+6*6 adds = 216 = 6*s^2 + 2s + 2 = 230
    for (i = 0; i < NWORDS; i++)
    {
        for (j = 0; j < i; j++)
        {
            // executed 0+1+2+3=6 times for nwords=4
            // executed 0+1+2+3+4=10 times for nwords=5
            // executed 0+1+2+3+4+5=15 times for nwords=6
            // each with 7 reads, 3 writes, 6 adds, 2 muls
            spMul2Acc(u->data[j], v->data[i - j], mdata->n->data[i - j], s->data[j], t);
        }
        // 5 reads, 3 writes, 3 adds, 1 mul
        spMulAddc(u->data[i], v->data[0], t);
        // 2 reads, 2 writes, 1 mul
        spMul(mdata->rho, t[0], &s->data[i], &c);
        // 5 reads, 3 writes, 3 adds, 1 mul
        spMulAddcr(mdata->n->data[0], s->data[i], t);
    }

    for (i = NWORDS; i < 2 * NWORDS; i++)
    {
        for (j = i - NWORDS + 1; j < NWORDS; j++)
        {
            // executed 0+1+2+3=6 times for nwords=4
            // each with 7 reads, 3 writes, 6 adds, 2 muls
            spMul2Acc(u->data[j], v->data[i - j], mdata->n->data[i - j], s->data[j], t);
        }
        // 3 reads, 3 writes
        s->data[i - NWORDS] = t[0];
        t[0] = t[1];
        t[1] = t[2];
        t[2] = 0;
    }

    // only get the result < 2^r, not strictly less than N
    if (t[0])
    {
        //printf("final subtraction: \n");
        mpSub(s->data, mdata->n->data, w->data, NWORDS);
    }
    else
    {
        memcpy(w->data, s->data, NWORDS * sizeof(uint64_t));
    }
    w->size = NWORDS;

    return;
}

void mulmod_bps(monty* mdata, bignum* u, bignum* v, bignum* w, bignum* s)
{
    // integrate multiply and reduction steps, alternating
    // between iterations of the inner loops while keeping
    // computations in product scanning form.
    int i, j, p;
    uint64_t acc[3] = { 0, 0, 0 };

    if (!initialized)
    {
        col = (uint64_t * *)malloc(BLOCKSIZE * sizeof(uint64_t*));
        for (i = 0; i < BLOCKSIZE; i++)
        {
            col[i] = (uint64_t*)malloc(3 * sizeof(uint64_t));
        }
        initialized = 1;
    }

    // compute total operations in 4 key categories and compare to FIPS 
    // figures from Koc's paper.  Muls are the same.  Slightly more adds
    // that are amortized towards Koc's figure as operand size increases.
    // Extra adds because the block accumulators must be added into the 
    // column accumulator for each block-column.  Reads and writes both
    // start out higher than FIPS but quickly amortize well below FIPS
    // as more and more full blocks are utilized.  Note that the reads/writes
    // could improve further with specialized routines for the partial blocks
    // that hold the block accumulator words in registers when possible.

    // total muls:
    // (nb*(nb - 1) * 2*B^2) + (nb * B*(B + 1)) + nb * B + (nb * B(B - 1)) muls
    // check nwords = 4: nb=1, B=4
    // 0 + 1*4(5) + 1*4(3) + 1*4 = 36 muls = 2*s^2 + s = 2*16 + 4 = 36 (100%)
    // check nwords = 8: nb=2, B=4
    // 2*1*2*4^2 + 2*4(5) + 2*4(3) + 2*4 = 136 muls = 2*s^2 + s = 2*64+8 = 136 (100%)
    // check nwords = 12: nb=3, B=4
    // 3*2*2*4^2 + 3*4(5) + 3*4(3) + 3*4 = 300 muls = 2*s^2 + s = 2*64+8 = 300 (100%)

    // total adds
    // (nb*(nb - 1) * 6*B^2) + (nb * 3*B*(B + 1)) + 2*nb*3*B + (nb * 3 * B(B - 1)) adds
    // check nwords = 4: nb=1, B=4
    // 0 + 1*3*4(5) + 1*3*4(3) + 2*1*3*4 = 120 adds = 6*s^2 + 2s + 2 = 6*16 + 24 = 106 (113%)
    // check nwords = 8: nb=2, B=4
    // 2*1*6*4^2 + 2*3*4(5) + 2*3*4(3) + 2*2*3*4 = 432 adds = 6*s^2 + 2s + 2 = 402 (107%)
    // check nwords = 12: nb=3, B=4
    // 3*2*6*4^2 + 3*3*4(5) + 3*3*4(3) + 3*2*3*4 = 936 adds = 6*s^2 + 2s + 2 = 890 (105%)

    // total reads
    // (nb*(nb - 1) * (4*B^2 + 3*B)) + (nb * 5*B*(B + 1)) + 2*nb*6*B + (nb * 5 * B(B - 1)) reads
    // check nwords = 4: nb=1, B=4
    // 0 + 1*5*4(5) + 1*5*4(3) + 2*1*6*4 = 208 reads; 9*s^2 + 8s + 2 = 9*16 + 32 + 2 = 178 (117%)
    // check nwords = 8: nb=2, B=4
    // 2*1*(4*4^2+3*4) + 2*5*4(5) + 2*5*4(3) + 2*2*6*4 = 568 reads = 9*8^2 + 8*8 + 2 = 642 (88%)
    // check nwords = 12: nb=3, B=4
    // 3*2*(4*4^2+3*4) + 3*5*4(5) + 3*5*4(3) + 3*2*6*4 = 1080 reads = 9*8^2 + 8*8 + 2 = 1394 (77%)

    // total writes
    // (nb*(nb - 1) * 3*B) + (nb * 3*B*(B + 1)) + 2*nb*4*B + (nb * 3 * B(B - 1)) writes
    // check nwords = 4: nb=1, B=4
    // 0 + 1*3*4(5) + 1*3*4(3) + 2*1*4*4 = 128 writes; 5*s^2 + 8s + 1 = 5*16 + 32 + 1 = 113 (113%)
    // check nwords = 8: nb=2, B=4
    // 2*1*(3*4) + 2*3*4(5) + 2*3*4(3) + 2*2*4*4 = 280 writes = 5*8^2 + 8*8 + 1 = 385 (73%)
    // check nwords = 12: nb=3, B=4
    // 3*2*(3*4) + 3*3*4(5) + 3*3*4(3) + 3*2*4*4 = 456 writes = 5*12^2 + 8*12 + 1 = 817(56%)

    for (i = 0; i < NBLOCKS; i++)
    {
        // 4 writes
        for (j = 0; j < BLOCKSIZE; j++)
        {
            col[j][0] = col[j][1] = col[j][2] = 0;
        }

        // nb*(nb - 1)/2 * B *(3 + 4*B reads, 3 writes, 6 * B adds, 2 * B muls)
        for (j = i; j > 0; j--)
        {
            // 0+1+2+3 iterations for nblocks=4
            // each with B*(3 + 4*B reads, 3 writes, 6 * B adds, 2 * B muls) = 76r, 12w, 96a, 32m
            //bps_fullblock_bsz4x2(u->data, v->data, s->data, mdata->n->data, i, j, col);
            bps_fullblock_bsz4(u->data, v->data, i, j, col);
            bps_fullblock_bsz4(s->data, mdata->n->data, i, j, col);
        }
        
        // nb * B*(B + 1)/2 * (5 reads, 3 writes, 3 adds, 1 mul)
        bps_pblock1(u->data, v->data, i, col);

        // nb * B*(B + 1)/2 * (5 reads, 3 writes, 3 adds, 1 mul)
        // plus nb * B * (6 reads, 3 adds, 1 mul, 4 writes)
        bps_final(s->data, mdata->n->data, col, i, acc, mdata->rho);
    }

    for (i = NBLOCKS; i < 2 * NBLOCKS; i++)
    {
        // 4 writes
        for (j = 0; j < BLOCKSIZE; j++)
        {
            col[j][0] = col[j][1] = col[j][2] = 0;
        }

        // nb*(nb - 1)/2 * B *(3 + 4*B reads, 3 writes, 6 * B adds, 2 * B muls)
        for (j = i - NBLOCKS + 1; j < NBLOCKS; j++)
        {
            //bps_fullblock_bsz4x2(u->data, v->data, s->data, mdata->n->data, i, j, col);
            bps_fullblock_bsz4(u->data, v->data, i, j, col);
            bps_fullblock_bsz4(s->data, mdata->n->data, i, j, col);
        }

        // nb * B(B - 1)/2 * 5 reads, 3 writes, 3 adds, 1 mul
        bps_pblock2(u->data, v->data, i, col);
        // nb * B(B - 1)/2 * 5 reads, 3 writes, 3 adds, 1 mul
        bps_pblock2(s->data, mdata->n->data, i, col);

        // nb * B* (3 adds + 6 reads, 4 writes)
        for (p = 0; p < BLOCKSIZE; p++)
        {
            int c = 0;
            c = _addcarry_u64(c, acc[0], col[p][0], &acc[0]);
            c = _addcarry_u64(c, acc[1], col[p][1], &acc[1]);
            c = _addcarry_u64(c, acc[2], col[p][2], &acc[2]);

            s->data[(i - NBLOCKS) * BLOCKSIZE + p] = acc[0];
            acc[0] = acc[1];
            acc[1] = acc[2];
            acc[2] = c;
        }
    }
    s->data[(i - NBLOCKS) * BLOCKSIZE] = 0;

    if (acc[0])
    {
        mpSub(s->data, mdata->n->data, w->data, NWORDS);
    }
    else
    {
        memcpy(w->data, s->data, NWORDS * sizeof(uint64_t));
    }
    w->size = NWORDS;

    return;
}


/* ==============================================================*/
/* Sqrmod routines that just call mulmod */
/* ==============================================================*/


void sqrmod_bps_mul(monty* mdata, bignum* u, bignum* w, bignum* s)
{
    mulmod_bps(mdata, u, u, w, s);
    return;
}

void sqrmod_sos_mul(monty *mdata, bignum * u, bignum * w, bignum * s)
{
    return mulmod_sos(mdata, u, u, w, s);
}

void sqrmod_fips_mul(monty *mdata, bignum * u, bignum * w, bignum * s)
{
    return mulmod_fips(mdata, u, u, w, s);
}

void sqrmod_fios_mul(monty *mdata, bignum * u, bignum * w, bignum * s)
{
    return mulmod_fios(mdata, u, u, w, s);
}

void sqrmod_cios_mul(monty *mdata, bignum * u, bignum * w, bignum * s)
{
    return mulmod_cios(mdata, u, u, w, s);
}


/* ==============================================================*/
/* Sqrmod routines */
/* ==============================================================*/


void sqrmod_sos(monty *mdata, bignum * u, bignum * w, bignum * s)
{
    // separate operand scanning. 
    // first square, then reduce
    int i, j;
    uint64_t c, p;

    memset(s->data, 0, (2 * NWORDS + 1) * sizeof(uint64_t));

    for (i = 0; i < NWORDS; i++)
    {
        spMulAdd2(u->data[i], u->data[i], s->data[2 * i], 0, &s->data[2 * i], &c);
        mpAdd1b(s->data + 2 * i + 1, c, s->data + 2 * i + 1, 2 * NWORDS);
        c = 0;

        for (j = i + 1; j < NWORDS; j++)
        {
            spMulDblAdd_1(u->data[j], u->data[i], c, &s->data[i + j], &c);
        }
        mpAdd1b(s->data + i + j + 1, c, s->data + i + j + 1, 2 * NWORDS - i - j);
    }
    s->size = 2 * NWORDS;

    for (i = 0; i < NWORDS; i++)
    {
        spMul(s->data[i], mdata->rho, &p, &c);
        c = 0;
        for (j = 0; j < NWORDS; j++)
        {
            spMulAdd2(mdata->n->data[j], p, s->data[i + j], c, &s->data[i + j], &c);
        }
        mpAdd1b(s->data + i + j, c, s->data + i + j, 2 * NWORDS - i - j);
    }
    s->size = 2 * NWORDS;

    // only get the result < 2^r, not strictly less than N
    if (s->data[2 * NWORDS])
    {
        mpSub(s->data + NWORDS, mdata->n->data, w->data, NWORDS);
    }
    else
    {
        memcpy(w->data, s->data + NWORDS, NWORDS * sizeof(uint64_t));
    }
    w->size = NWORDS;

    return;
}

void sqrmod_fips(monty *mdata, bignum * u, bignum * w, bignum * s)
{
    // integrate square and reduction steps, alternating
    // between iterations of the inner loops while keeping
    // computations in product scanning form.
    int i, j;
    uint64_t c, p;
    uint64_t t[3] = { 0, 0, 0 };

    memset(s->data, 0, (NWORDS + 1) * sizeof(uint64_t));

    for (i = 0; i < NWORDS; i++)
    {
        // terms that are doubled
        for (j = 0; j < i / 2; j++)
        {
            spSqrMulAcc(u->data[j], u->data[i - j], mdata->n->data[i - j], s->data[j], t);
        }

        if (i & 1)
        {
            spMulDblAdd_3(u->data[j], u->data[i - j], t);
        }
        else
        {
            spMulAddc(u->data[j], u->data[j], t);
        }

        // finish contribution from n (not a squaring operation)
        for (; j < i; j++)
        {
            spMulAddc(mdata->n->data[i - j], s->data[j], t);
        }

        // form this 'm' term
        spMul(mdata->rho, t[0], &s->data[i], &c);

        // final 't' update using this new 'm'
        spMulAddcr(mdata->n->data[0], s->data[i], t);
    }

    for (i = NWORDS; i < 2 * NWORDS - 1; i++)
    {
        for (j = i - NWORDS + 1; j < i / 2; j++)
        {
            spSqrMulAcc(u->data[j], u->data[i - j], mdata->n->data[i - j], s->data[j], t);
        }

        if (i & 1)
        {
            // no square term
            spMulDblAdd_3(u->data[j], u->data[i - j], t);
        }
        else
        {
            // square term
            spMulAddc(u->data[j], u->data[j], t);
        }

        // finish contribution from n (not a squaring operation)
        for (; j < NWORDS; j++)
        {
            spMulAddc(mdata->n->data[i - j], s->data[j], t);
        }

        s->data[i - NWORDS] = t[0];
        t[0] = t[1];
        t[1] = t[2];
        t[2] = 0;
    }
    s->data[i - NWORDS] = t[0];

    // only get the result < 2^r, not strictly less than N
    if (t[1])
    {
        mpSub(s->data, mdata->n->data, w->data, NWORDS);
    }
    else
    {
        memcpy(w->data, s->data, NWORDS * sizeof(uint64_t));
    }
    w->size = NWORDS;

    return;
}

void sqrmod_fios(monty *mdata, bignum * u, bignum * w, bignum * s)
{
    // integrate multiply and reduction steps, alternating
    // between iterations of the inner loops.
    int i, j;
    uint64_t c, c2, c3, p;

    memset(s->data, 0, (NWORDS + 1) * sizeof(uint64_t));

    i = 0;
    {
        uint64_t s_in = s->data[0];
        spMulAdd1(u->data[0], u->data[i], s->data[0], &s->data[0], &c);
        spAdd(s->data[1], c, s->data + 1, &c2);     // c2 applies to j+2 (next j+1)
        spMul(s->data[0], mdata->rho, &p, &c);
        spMulAdd1(mdata->n->data[0], p, s->data[0], &s->data[0], &c);

        // double all of the i=0 cross-terms.  then can eliminate them from 
        // subsequent i.
        s->data[NWORDS + 1] = 0;
        c3 = 0;
        for (j = 1; j < NWORDS; j++)
        {
            spMulAdd2(u->data[j], u->data[i], s->data[j], c, &s->data[j - 1], &c);
            spAdd3(s->data[j + 1], c, c2, s->data + j + 1, &c2);     // c2 applies to j+2 (next j+1)
            spMulAdd2(u->data[j], u->data[i], s->data[j - 1], 0, &s->data[j - 1], &c);
            spAdd3(s->data[j + 1], c, c3, s->data + j + 1, &c3);     // c3 applies to j+2 (next j+1)
            spMulAdd1(mdata->n->data[j], p, s->data[j - 1], &s->data[j - 1], &c);
        }
        spAdd(s->data[j], c, &s->data[j - 1], &s->data[j]);
        // the doubling of cross-terms means that sometimes we'll get 
        // an extra carry that needs to be accumulated.
        s->data[j] += (s->data[j + 1] + c3);
    }

    for (i = 1; i < NWORDS; i++)
    {
        // the i=0 case included all of these terms already...
        uint64_t s_in = s->data[0];
        spMul(s->data[0], mdata->rho, &p, &c);
        spMulAdd1(mdata->n->data[0], p, s->data[0], &s->data[0], &c);

        // skip to the square case.  previous outer loop iterations have already
        // included prior terms.
        for (j = 1; j < i; j++)
        {
            spMulAdd2(mdata->n->data[j], p, s->data[j], c, &s->data[j - 1], &c);
        }

        // square case
        spMulAdd2(u->data[j], u->data[i], s->data[j], c, &s->data[j - 1], &c);
        spAdd(s->data[j + 1], c, s->data + j + 1, &c2);     // c2 applies to j+2 (next j+1)
        spMulAdd1(mdata->n->data[j], p, s->data[j - 1], &s->data[j - 1], &c);
        j++;

        // double remaining terms
        s->data[NWORDS + 1] = 0;
        c3 = 0;
        for (; j < NWORDS; j++)
        {
            uint64_t t0, t1;
            spMul(u->data[j], u->data[i], &t0, &t1);
            spAdd3(t0, t0, c, &t0, &c);                     // c  applies to j+1 (next j)
            spMulAdd2(mdata->n->data[j], p, s->data[j], t0, &s->data[j - 1], &c3);   // c3 applies to j+1 (next j)            
            spAdd3(t1, t1, c2, &t1, &c2);                   // c2 applies to j+2 (next j+1)    
            spAdd3(s->data[j + 1], t1, c3, &s->data[j + 1], &c3);
            c2 += c3;
        }
        spAdd(s->data[j], c, &s->data[j - 1], &s->data[j]);
        // the doubling of cross-terms means that sometimes we'll get 
        // an extra carry that needs to be accumulated.
        s->data[j] += (s->data[j + 1] + c2);
    }

    // only get the result < 2^r, not strictly less than N
    if (s->data[NWORDS])
    {
        mpSub(s->data, mdata->n->data, w->data, NWORDS);
    }
    else
    {
        memcpy(w->data, s->data, NWORDS * sizeof(uint64_t));
    }
    w->size = NWORDS;

    return;
}

void sqrmod_cios(monty *mdata, bignum * u, bignum * w, bignum * s)
{
    // integrate multiply and reduction steps, alternating
    // between iterations of the outer loops.
    int i, j;
    uint64_t c, p;

    memset(s->data, 0, (NWORDS + 1) * sizeof(uint64_t));

    for (i = 0; i < NWORDS; i++)
    {
        // jump to square term.  skipped terms have already been included
        // in the result by previous doubled iterations.
        j = i;
        spMulAdd1(u->data[j], u->data[i], s->data[j], &s->data[j], &c);
        mpAdd1b(s->data + i + 1, c, s->data + i + 1, NWORDS - i - 1);

        // double the rest of the cross terms
        for (j++; j < NWORDS; j++)
        {
            spMulDblAdd_3(u->data[j], u->data[i], &s->data[j]);
        }

        spMul(s->data[0], mdata->rho, &p, &c);

        spMulAdd1(mdata->n->data[0], p, s->data[0], &s->data[0], &c);
        for (j = 1; j < NWORDS; j++)
        {
            spMulAdd2x(mdata->n->data[j], p, s->data[j], c, &s->data[j - 1], &c);
        }
        spAdd(s->data[j], c, &s->data[j - 1], &s->data[j]);
        // possible extra carry due to doubling of cross-terms.
        s->data[j] += s->data[j + 1];
        s->data[j + 1] = 0;
    }

    s->size = NWORDS;

    // only get the result < 2^r, not strictly less than N
    if (s->data[NWORDS])
    {
        mpSub(s->data, mdata->n->data, w->data, NWORDS);
    }
    else
    {
        memcpy(w->data, s->data, NWORDS * sizeof(uint64_t));
    }
    w->size = NWORDS;

    return;
}

void sqrmod_bps(monty* mdata, bignum* u, bignum* w, bignum* s)
{
    // integrate multiply and reduction steps, alternating
    // between iterations of the inner loops while keeping
    // computations in product scanning form.
    int i, j, p;
    uint64_t acc[3] = { 0, 0, 0 };
    bignum* v = u;

    if (!initialized)
    {

        col = (uint64_t * *)malloc(BLOCKSIZE * sizeof(uint64_t*));
        for (i = 0; i < BLOCKSIZE; i++)
        {
            col[i] = (uint64_t*)malloc(3 * sizeof(uint64_t));
        }
        initialized = 1;
    }

    for (i = 0; i < NBLOCKS; i++)
    {
        for (j = 0; j < BLOCKSIZE; j++)
        {
            col[j][0] = col[j][1] = col[j][2] = 0;
        }

        for (j = i; j > (i + 1) / 2; j--)
        {
            bps_fullblock_bsz4x4(u->data, v->data, i, j, col);
        }

        if (i & 1)
        {
            bps_pblock1p_bsz4(u->data, v->data, i, col);
        }
        else
        {
            bps_pblock1pp_bsz4(u->data, v->data, i, col);
        }

        for (j = 0; j < BLOCKSIZE; j++)
        {
            // double each accumulator now that all u * v terms are accumulated
            uint64_t hibit = col[j][0] >> 63;
            uint64_t hibit2 = col[j][1] >> 63;
            col[j][0] <<= 1;
            col[j][1] = (col[j][1] << 1) + hibit;
            col[j][2] = (col[j][2] << 1) + hibit2;
        }

        // accumulate the square terms (not doubled)
        spMulAddc(u->data[i * 2 + 0], u->data[i * 2 + 0], col[0]);
        spMulAddc(u->data[i * 2 + 1], u->data[i * 2 + 1], col[2]);

        // accumulate the s * n terms (not doubled)
        for (j = i; j > 0; j--)
        {
            bps_fullblock_bsz4x4(s->data, mdata->n->data, i, j, col);
        }
        bps_final(s->data, mdata->n->data, col, i, acc, mdata->rho);
    }

    for (i = NBLOCKS; i < 2 * NBLOCKS; i++)
    {
        for (j = 0; j < BLOCKSIZE; j++)
        {
            col[j][0] = col[j][1] = col[j][2] = 0;
        }

        for (j = i - NBLOCKS + 1; j < NBLOCKS; j++)
        {
            // Compare the max index of B to the min index of A 
            // for this potential block
            if (((j - 1) * BLOCKSIZE + (2 * BLOCKSIZE - 1)) > ((i - j) * BLOCKSIZE))
            {
                // if max(Bid) > min(Aid) then this block has doubles: it is not a full block
                break;
            }
            bps_fullblock_bsz4x4(u->data, v->data, i, j, col);
        }

        if (i & 1)
        {
            bps_pblock2p_bsz4(u->data, v->data, i, col);
        }
        else
        {
            bps_pblock2pp_bsz4(u->data, v->data, i, col);
        }

        for (j = 0; j < BLOCKSIZE; j++)
        {
            // double each accumulator now that all u * v terms are accumulated
            uint64_t hibit = col[j][0] >> 63;
            uint64_t hibit2 = col[j][1] >> 63;
            col[j][0] <<= 1;
            col[j][1] = (col[j][1] << 1) + hibit;
            col[j][2] = (col[j][2] << 1) + hibit2;
        }

        // accumulate the square terms (not doubled)
        spMulAddc(u->data[i * 2 + 0], u->data[i * 2 + 0], col[0]);
        spMulAddc(u->data[i * 2 + 1], u->data[i * 2 + 1], col[2]);

        // accumulate the s * n terms (not doubled)
        for (j = i - NBLOCKS + 1; j < NBLOCKS; j++)
        {
            bps_fullblock_bsz4x4(s->data, mdata->n->data, i, j, col);

        }
        bps_pblock2(s->data, mdata->n->data, i, col);

        for (p = 0; p < BLOCKSIZE; p++)
        {
            int c = 0;
            c = _addcarry_u64(c, acc[0], col[p][0], &acc[0]);
            c = _addcarry_u64(c, acc[1], col[p][1], &acc[1]);
            _addcarry_u64(c, acc[2], col[p][2], &acc[2]);

            s->data[(i - NBLOCKS) * BLOCKSIZE + p] = acc[0];
            acc[0] = acc[1];
            acc[1] = acc[2];
            acc[2] = 0;
        }
    }

    if (acc[0])
    {
        mpSub(s->data, mdata->n->data, w->data, NWORDS);
    }
    else
    {
        memcpy(w->data, s->data, NWORDS * sizeof(uint64_t));
    }
    w->size = NWORDS;

    return;
}
