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
#include "gmp.h"
#include "omp.h"

uint64_t *LCG_STATE;

uint64_t spRand64(uint64_t *state)
{
    // advance the state of the LCG and return the appropriate result.
    // assume lower = 0 and upper = maxint
    *state = 6364136223846793005ULL * (*state) + 1442695040888963407ULL;
    return *state;
}

// FNV-1 hash algorithm:
// http://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
uint64_t hash64(uint64_t in)
{
    uint64_t hash = 14695981039346656037ULL;
    uint64_t prime = 1099511628211ULL;
    uint64_t hash_mask;
    uint64_t xor;

    hash = hash * prime;
    hash_mask = 0xffffffffffffff00ULL;
    xor = hash ^ in;
    hash = (hash & hash_mask) | (xor & (~hash_mask));

    hash = hash * prime;
    hash_mask = 0xffffffffffff00ffULL;
    xor = hash ^ in;
    hash = (hash & hash_mask) | (xor & (~hash_mask));

    hash = hash * prime;
    hash_mask = 0xffffffffff00ffffULL;
    xor = hash ^ in;
    hash = (hash & hash_mask) | (xor & (~hash_mask));

    hash = hash * prime;
    hash_mask = 0xffffffff00ffffffULL;
    xor = hash ^ in;
    hash = (hash & hash_mask) | (xor & (~hash_mask));

    hash = hash * prime;
    hash_mask = 0xffffff00ffffffffULL;
    xor = hash ^ in;
    hash = (hash & hash_mask) | (xor & (~hash_mask));

    hash = hash * prime;
    hash_mask = 0xffff00ffffffffffULL;
    xor = hash ^ in;
    hash = (hash & hash_mask) | (xor & (~hash_mask));

    hash = hash * prime;
    hash_mask = 0xff00ffffffffffffULL;
    xor = hash ^ in;
    hash = (hash & hash_mask) | (xor & (~hash_mask));

    hash = hash * prime;
    hash_mask = 0x00ffffffffffffffULL;
    xor = hash ^ in;
    hash = (hash & hash_mask) | (xor & (~hash_mask));

    return hash;
}

double my_difftime(struct timeval * start, struct timeval * end)
{
    double secs;
    double usecs;

    if (start->tv_sec == end->tv_sec) {
        secs = 0;
        usecs = end->tv_usec - start->tv_usec;
    }
    else {
        usecs = 1000000 - start->tv_usec;
        secs = end->tv_sec - (start->tv_sec + 1);
        usecs += end->tv_usec;
        if (usecs >= 1000000) {
            usecs -= 1000000;
            secs += 1;
        }
    }

    return secs + usecs / 1000000.;
}

void extract_bignum_from_vec_to_mpz(mpz_t dest, bignum *vec_src, int num, int sz)
{
    int j;

    if (dest == NULL)
    {
        printf("invalid dest address in extract_vec_bignum_from_vec_to_mpz\n");
    }

    mpz_set_ui(dest, 0);
    for (j = sz - 1; j >= 0; j--)
    {
        mpz_mul_2exp(dest, dest, DIGITBITS);
        mpz_add_ui(dest, dest, vec_src->data[num + j * VECLEN]);
    }

    return;
}

void vecpmodtest(int do_verification, int threads, int verbose)
{
    // test the pmod by comparing all results to those computed using
    // validated scalar code.
    double *elapsed_time;
    int t;
    //gmp_randstate_t rng_state;

    //gmp_randinit_default(rng_state);
    elapsed_time = (double *)malloc(threads * sizeof(double));

    LCG_STATE = (uint64_t *)malloc(threads * sizeof(uint64_t));

    for (t = 0; t < threads; t++)
    {
        LCG_STATE[t] = hash64(t);
    }

    printf("commencing test: all variable (random)\n");
#pragma omp parallel num_threads(threads)
    {
        int i, j;

        // timing variables
        struct timeval stopt;	// stop time of this job
        struct timeval startt;	// start time of this job
        double t_time;

        mpz_t base, exp, mod, t1, t2;
        
        int loc_iterations;
        int tid = omp_get_thread_num();
        monty *mtest;

        // vector bignums
        bignum *b = vecInit();
        bignum *d = vecInit();
        bignum *m = vecInit();
        bignum *e = vecInit();
        bignum *s = vecInit();
        bignum *one = vecInit();

        mpz_init(base);
        mpz_init(exp);
        mpz_init(mod);
        mpz_init(t1);
        mpz_init(t2);
        
        //gmp_randseed_ui(rng_state, tid);

        // attempt to scale the number of iterations with input size
        // so this doesn't take forever.
        loc_iterations = 100000 * 2 / (NWORDS * DIGITBITS);

        if (MAXBITS >= 4096)
            loc_iterations *= 1; 
        else if (MAXBITS >= 2048)
            loc_iterations *= 2;
        else if (MAXBITS >= 1024)
            loc_iterations *= 5;
        else if (MAXBITS >= 512)
            loc_iterations *= 10;
        else if (MAXBITS >= 256)
            loc_iterations *= 25;
        else
            loc_iterations *= 100;

#ifdef BASE52
        loc_iterations *= 3;
#endif

#ifdef TARGET_KNL
        loc_iterations /= 3;
#endif

        mtest = monty_alloc();

#pragma omp barrier

        gettimeofday(&startt, NULL);
        
        for (j = 0; j < VECLEN; j++)
        {
            one->data[j] = 1;
        }

        printf("thread %d starting %d iterations\n", tid, loc_iterations);

        // now do the calculation "b^e % m" a bunch of times
        for (i = 0; i < loc_iterations; i++)
        {

#ifdef BASE52
            //int tmp = ceil(MAXBITS / 64);
            memset(m->data, 0, MAXBITS * 2 * VECLEN / 8);
            memset(e->data, 0, MAXBITS * 2 * VECLEN / 8);
            memset(b->data, 0, MAXBITS * 2 * VECLEN / 8);

            for (j = 0; j < VECLEN; j++)
            {
                int k;
                for (k = 0; k < NWORDS; k++)
                {
                    uint64_t r1 = spRand64(&LCG_STATE[t]);
                    uint64_t r2 = spRand64(&LCG_STATE[t]);
                    uint64_t r3 = spRand64(&LCG_STATE[t]);

                    m->data[k * VECLEN + j] = r1 & MAXDIGIT;
                    b->data[k * VECLEN + j] = r2 & MAXDIGIT;
                    e->data[k * VECLEN + j] = r3 & MAXDIGIT;
                }
            }

#else
            memset(m->data, 0, MAXBITS * 2 * VECLEN / 8);
            memset(e->data, 0, MAXBITS * 2 * VECLEN / 8);
            memset(b->data, 0, MAXBITS * 2 * VECLEN / 8);

            for (j = 0; j < VECLEN; j++)
            {
                int k;
                for (k = 0; k < NWORDS; k++)
                {
                    uint64_t r1 = spRand64(&LCG_STATE[t]);
                    uint64_t r2 = spRand64(&LCG_STATE[t]);
                    uint64_t r3 = spRand64(&LCG_STATE[t]);

                    m->data[k * VECLEN + j] = r1 & MAXDIGIT;
                    b->data[k * VECLEN + j] = r2 & MAXDIGIT;
                    e->data[k * VECLEN + j] = r3 & MAXDIGIT;
                }
            }
#endif
            for (j = 0; j < VECLEN; j++)
                m->data[j] |= 0x1;

            if (0)
            {
                continue;
            }

            if (verbose > 1)
            {
                for (j = 0; j < VECLEN; j++)
                {
                    extract_bignum_from_vec_to_mpz(base, b, j, NWORDS);
                    extract_bignum_from_vec_to_mpz(exp, e, j, NWORDS);
                    extract_bignum_from_vec_to_mpz(mod, m, j, NWORDS);

                    gmp_printf("init%d:\n\tbase = %Zx\n\texp = %Zx\n\tmod = %Zx\n",
                        j, base, exp, mod);
                }
            }

            // now we actually do the (vectorized) montgomery initialization
            // on our vector of random moduli.
            monty_init_vec(mtest, m, 0);

            if (verbose > 1)
            {
                for (j = 0; j < VECLEN; j++)
                {
                    extract_bignum_from_vec_to_mpz(base, mtest->rhat, j, NWORDS);
                    extract_bignum_from_vec_to_mpz(mod, mtest->one, j, NWORDS);

                    gmp_printf("init%d:\n\trhat = %Zx\n\tone = %Zx\n\trho = %08x\n",
                        j, base, mod, mtest->vrho[j]);
                }
            }

            vecmulmod_ptr(b, mtest->rhat, b, m, s, mtest);      // monty rep

            if (verbose > 1)
            {
                for (j = 0; j < VECLEN; j++)
                {
                    extract_bignum_from_vec_to_mpz(base, b, j, NWORDS);

                    gmp_printf("monty(base%d) = %Zx\n", j, base);
                }
            }

            vecmodexp_ptr(d, b, e, m, s, mtest->one, mtest);    // powm
            vecmulmod_ptr(d, one, d, m, s, mtest);              // normal rep

            if (verbose > 1)
            {
                for (j = 0; j < VECLEN; j++)
                {
                    extract_bignum_from_vec_to_mpz(base, d, j, NWORDS);

                    gmp_printf("modexp%d = %Zx\n", j, base);
                }
            }

            // now verify each result
            if (do_verification)
            {
                vecmulmod_ptr(b, one, b, m, s, mtest);              // normal rep
                for (j = 0; j < VECLEN; j++)
                {
                    extract_bignum_from_vec_to_mpz(t1, d, j, NWORDS);
                    extract_bignum_from_vec_to_mpz(base, b, j, NWORDS);
                    extract_bignum_from_vec_to_mpz(exp, e, j, NWORDS);
                    extract_bignum_from_vec_to_mpz(mod, m, j, NWORDS);

                    mpz_powm(t2, base, exp, mod);

                    if (verbose)
                    {
                        gmp_printf("iteration %d lane %d:\n\tgmp = %Zx\n\ttest = %Zx\n",
                            i, j, t2, t1);
                    }

                    if (mpz_cmp(t1, t2) != 0)
                    {
                        gmp_printf("iteration %d error lane %d:\nbase = %Zx\nexp = %Zx\nmod = %Zx\ngmp = %Zx\ntest = %Zx\n",
                            i, j, base, exp, mod, t2, t1);
                        exit(1);
                    }

                }
            }           
        }

        monty_free(mtest);

        if ((tid == 0) && (do_verification == 1))
            printf("verified %d x 16 vecModExp (all variable) results\n", loc_iterations);

        gettimeofday(&stopt, NULL);
        t_time = my_difftime(&startt, &stopt);
        elapsed_time[tid] = t_time;

        if (tid == 0)
            printf("Test with %d iterations took %1.4f seconds.\n", loc_iterations, t_time);

        mpz_clear(t1);
        mpz_clear(t2);
        mpz_clear(base);
        mpz_clear(mod);
        mpz_clear(exp);
        vecFree(m);
        vecFree(b);
        vecFree(d);
        vecFree(s);
        vecFree(e);
    }

    {
        int i;
        double sum = 0.0;
        double min_t = 9999999999.;
        double max_t = 0.;

        for (i = 0; i < threads; i++)
        {
            sum += elapsed_time[i];
            if (elapsed_time[i] < min_t)
                min_t = elapsed_time[i];
            if (elapsed_time[i] > max_t)
                max_t = elapsed_time[i];
        }

        printf("average elapsed time = %1.4f\n", sum / threads);
        printf("min elapsed time = %1.4f\n", min_t);
        printf("max elapsed time = %1.4f\n", max_t);
    }

    free(elapsed_time);
    free(LCG_STATE);

    printf("\n\n");

    return;
}

void vecmultest(int do_verification, int threads, int verbose)
{
    // test the pmod by comparing all results to those computed using
    // validated scalar code.
    double* elapsed_time;
    int t;
    //gmp_randstate_t rng_state;

    //gmp_randinit_default(rng_state);
    elapsed_time = (double*)malloc(threads * sizeof(double));

    LCG_STATE = (uint64_t*)malloc(threads * sizeof(uint64_t));

    for (t = 0; t < threads; t++)
    {
        LCG_STATE[t] = hash64(t);
    }

    do_verification = 0;
    printf("commencing test mulmod: all variable (random)\n");
#pragma omp parallel num_threads(threads)
    {
        int i, j;

        // timing variables
        struct timeval stopt;	// stop time of this job
        struct timeval startt;	// start time of this job
        double t_time;

        mpz_t base, exp, mod, t1, t2;

        int loc_iterations;
        int tid = omp_get_thread_num();
        monty* mtest;

        // vector bignums
        bignum* b = vecInit();
        bignum* d = vecInit();
        bignum* m = vecInit();
        bignum* e = vecInit();
        bignum* s = vecInit();
        bignum* one = vecInit();

        mpz_init(base);
        mpz_init(exp);
        mpz_init(mod);
        mpz_init(t1);
        mpz_init(t2);

        //gmp_randseed_ui(rng_state, tid);

        // attempt to scale the number of iterations with input size
        // so this doesn't take forever.
        loc_iterations = 100000 * 2 / (NWORDS * DIGITBITS);

        if (MAXBITS >= 4096)
            loc_iterations *= 1;
        else if (MAXBITS >= 2048)
            loc_iterations *= 2;
        else if (MAXBITS >= 1024)
            loc_iterations *= 5;
        else if (MAXBITS >= 512)
            loc_iterations *= 10;
        else if (MAXBITS >= 256)
            loc_iterations *= 25;
        else
            loc_iterations *= 100;

#ifdef BASE52
        loc_iterations *= 3;
#endif

#ifdef TARGET_KNL
        loc_iterations /= 3;
#endif

        loc_iterations *= 20000;
        mtest = monty_alloc();

#pragma omp barrier

        gettimeofday(&startt, NULL);

        for (j = 0; j < VECLEN; j++)
        {
            one->data[j] = 1;
        }

#ifdef BASE52
        //int tmp = ceil(MAXBITS / 64);
        memset(m->data, 0, MAXBITS * 2 * VECLEN / 8);
        memset(e->data, 0, MAXBITS * 2 * VECLEN / 8);
        memset(b->data, 0, MAXBITS * 2 * VECLEN / 8);

        for (j = 0; j < VECLEN; j++)
        {
            int k;
            for (k = 0; k < NWORDS; k++)
            {
                uint64_t r1 = spRand64(&LCG_STATE[t]);
                uint64_t r2 = spRand64(&LCG_STATE[t]);
                uint64_t r3 = spRand64(&LCG_STATE[t]);

                m->data[k * VECLEN + j] = r1 & MAXDIGIT;
                b->data[k * VECLEN + j] = r2 & MAXDIGIT;
                e->data[k * VECLEN + j] = r3 & MAXDIGIT;
            }
        }

#else
        memset(m->data, 0, MAXBITS * 2 * VECLEN / 8);
        memset(e->data, 0, MAXBITS * 2 * VECLEN / 8);
        memset(b->data, 0, MAXBITS * 2 * VECLEN / 8);

        for (j = 0; j < VECLEN; j++)
        {
            int k;
            for (k = 0; k < NWORDS; k++)
            {
                uint64_t r1 = spRand64(&LCG_STATE[t]);
                uint64_t r2 = spRand64(&LCG_STATE[t]);
                uint64_t r3 = spRand64(&LCG_STATE[t]);

                m->data[k * VECLEN + j] = r1 & MAXDIGIT;
                b->data[k * VECLEN + j] = r2 & MAXDIGIT;
                e->data[k * VECLEN + j] = r3 & MAXDIGIT;
            }
        }
#endif
        for (j = 0; j < VECLEN; j++)
            m->data[j] |= 0x1;

        if (verbose > 1)
        {
            for (j = 0; j < VECLEN; j++)
            {
                extract_bignum_from_vec_to_mpz(base, b, j, NWORDS);
                extract_bignum_from_vec_to_mpz(exp, e, j, NWORDS);
                extract_bignum_from_vec_to_mpz(mod, m, j, NWORDS);

                gmp_printf("init%d:\n\tbase = %Zx\n\texp = %Zx\n\tmod = %Zx\n",
                    j, base, exp, mod);
            }
        }

        // now we actually do the (vectorized) montgomery initialization
        // on our vector of random moduli.
        monty_init_vec(mtest, m, 0);

        if (verbose > 1)
        {
            for (j = 0; j < VECLEN; j++)
            {
                extract_bignum_from_vec_to_mpz(base, mtest->rhat, j, NWORDS);
                extract_bignum_from_vec_to_mpz(mod, mtest->one, j, NWORDS);

                gmp_printf("init%d:\n\trhat = %Zx\n\tone = %Zx\n\trho = %08x\n",
                    j, base, mod, mtest->vrho[j]);
            }
        }

        vecmulmod_ptr(b, mtest->rhat, b, m, s, mtest);      // monty rep
        vecmulmod_ptr(e, mtest->rhat, e, m, s, mtest);      // monty rep

        printf("thread %d starting %d iterations\n", tid, loc_iterations);

        // now do the calculation "b^e % m" a bunch of times
        for (i = 0; i < loc_iterations; i++)
        {
            if (verbose > 1)
            {
                for (j = 0; j < VECLEN; j++)
                {
                    extract_bignum_from_vec_to_mpz(base, b, j, NWORDS);

                    gmp_printf("monty(base%d) = %Zx\n", j, base);
                }
            }

            vecmulmod_ptr(b, e, b, m, s, mtest);

            if (verbose > 1)
            {
                for (j = 0; j < VECLEN; j++)
                {
                    extract_bignum_from_vec_to_mpz(base, d, j, NWORDS);

                    gmp_printf("modexp%d = %Zx\n", j, base);
                }
            }

            // now verify each result
            if (do_verification)
            {
                vecmulmod_ptr(b, one, b, m, s, mtest);              // normal rep
                for (j = 0; j < VECLEN; j++)
                {
                    extract_bignum_from_vec_to_mpz(t1, d, j, NWORDS);
                    extract_bignum_from_vec_to_mpz(base, b, j, NWORDS);
                    extract_bignum_from_vec_to_mpz(exp, e, j, NWORDS);
                    extract_bignum_from_vec_to_mpz(mod, m, j, NWORDS);

                    mpz_powm(t2, base, exp, mod);

                    if (verbose)
                    {
                        gmp_printf("iteration %d lane %d:\n\tgmp = %Zx\n\ttest = %Zx\n",
                            i, j, t2, t1);
                    }

                    if (mpz_cmp(t1, t2) != 0)
                    {
                        gmp_printf("iteration %d error lane %d:\nbase = %Zx\nexp = %Zx\nmod = %Zx\ngmp = %Zx\ntest = %Zx\n",
                            i, j, base, exp, mod, t2, t1);
                        exit(1);
                    }

                }
            }
        }

        monty_free(mtest);

        if ((tid == 0) && (do_verification == 1))
            printf("verified %d x 16 vecModExp (all variable) results\n", loc_iterations);

        gettimeofday(&stopt, NULL);
        t_time = my_difftime(&startt, &stopt);
        elapsed_time[tid] = t_time;

        if (tid == 0)
            printf("Test with %d iterations took %1.4f seconds.\n", loc_iterations, t_time);

        mpz_clear(t1);
        mpz_clear(t2);
        mpz_clear(base);
        mpz_clear(mod);
        mpz_clear(exp);
        vecFree(m);
        vecFree(b);
        vecFree(d);
        vecFree(s);
        vecFree(e);
    }

    {
        int i;
        double sum = 0.0;
        double min_t = 9999999999.;
        double max_t = 0.;

        for (i = 0; i < threads; i++)
        {
            sum += elapsed_time[i];
            if (elapsed_time[i] < min_t)
                min_t = elapsed_time[i];
            if (elapsed_time[i] > max_t)
                max_t = elapsed_time[i];
        }

        printf("average elapsed time = %1.4f\n", sum / threads);
        printf("min elapsed time = %1.4f\n", min_t);
        printf("max elapsed time = %1.4f\n", max_t);
    }

    free(elapsed_time);
    free(LCG_STATE);

    printf("\n\n");

    return;
}

void vecsqrtest(int do_verification, int threads, int verbose)
{
    // test the pmod by comparing all results to those computed using
    // validated scalar code.
    double* elapsed_time;
    int t;
    //gmp_randstate_t rng_state;

    //gmp_randinit_default(rng_state);
    elapsed_time = (double*)malloc(threads * sizeof(double));

    LCG_STATE = (uint64_t*)malloc(threads * sizeof(uint64_t));

    for (t = 0; t < threads; t++)
    {
        LCG_STATE[t] = hash64(t);
    }

    do_verification = 0;
    printf("commencing test sqrmod: all variable (random)\n");
#pragma omp parallel num_threads(threads)
    {
        int i, j;

        // timing variables
        struct timeval stopt;	// stop time of this job
        struct timeval startt;	// start time of this job
        double t_time;

        mpz_t base, exp, mod, t1, t2;

        int loc_iterations;
        int tid = omp_get_thread_num();
        monty* mtest;

        // vector bignums
        bignum* b = vecInit();
        bignum* d = vecInit();
        bignum* m = vecInit();
        bignum* e = vecInit();
        bignum* s = vecInit();
        bignum* one = vecInit();

        mpz_init(base);
        mpz_init(exp);
        mpz_init(mod);
        mpz_init(t1);
        mpz_init(t2);

        //gmp_randseed_ui(rng_state, tid);

        // attempt to scale the number of iterations with input size
        // so this doesn't take forever.
        loc_iterations = 100000 * 2 / (NWORDS * DIGITBITS);

        if (MAXBITS >= 4096)
            loc_iterations *= 1;
        else if (MAXBITS >= 2048)
            loc_iterations *= 2;
        else if (MAXBITS >= 1024)
            loc_iterations *= 5;
        else if (MAXBITS >= 512)
            loc_iterations *= 10;
        else if (MAXBITS >= 256)
            loc_iterations *= 25;
        else
            loc_iterations *= 100;

#ifdef BASE52
        loc_iterations *= 3;
#endif

#ifdef TARGET_KNL
        loc_iterations /= 3;
#endif

        loc_iterations *= 20000;
        mtest = monty_alloc();

#pragma omp barrier

        gettimeofday(&startt, NULL);

        for (j = 0; j < VECLEN; j++)
        {
            one->data[j] = 1;
        }

#ifdef BASE52
        //int tmp = ceil(MAXBITS / 64);
        memset(m->data, 0, MAXBITS * 2 * VECLEN / 8);
        memset(e->data, 0, MAXBITS * 2 * VECLEN / 8);
        memset(b->data, 0, MAXBITS * 2 * VECLEN / 8);

        for (j = 0; j < VECLEN; j++)
        {
            int k;
            for (k = 0; k < NWORDS; k++)
            {
                uint64_t r1 = spRand64(&LCG_STATE[t]);
                uint64_t r2 = spRand64(&LCG_STATE[t]);
                uint64_t r3 = spRand64(&LCG_STATE[t]);

                m->data[k * VECLEN + j] = r1 & MAXDIGIT;
                b->data[k * VECLEN + j] = r2 & MAXDIGIT;
                e->data[k * VECLEN + j] = r3 & MAXDIGIT;
            }
        }

#else
        memset(m->data, 0, MAXBITS * 2 * VECLEN / 8);
        memset(e->data, 0, MAXBITS * 2 * VECLEN / 8);
        memset(b->data, 0, MAXBITS * 2 * VECLEN / 8);

        for (j = 0; j < VECLEN; j++)
        {
            int k;
            for (k = 0; k < NWORDS; k++)
            {
                uint64_t r1 = spRand64(&LCG_STATE[t]);
                uint64_t r2 = spRand64(&LCG_STATE[t]);
                uint64_t r3 = spRand64(&LCG_STATE[t]);

                m->data[k * VECLEN + j] = r1 & MAXDIGIT;
                b->data[k * VECLEN + j] = r2 & MAXDIGIT;
                e->data[k * VECLEN + j] = r3 & MAXDIGIT;
            }
        }
#endif
        for (j = 0; j < VECLEN; j++)
            m->data[j] |= 0x1;

        if (verbose > 1)
        {
            for (j = 0; j < VECLEN; j++)
            {
                extract_bignum_from_vec_to_mpz(base, b, j, NWORDS);
                extract_bignum_from_vec_to_mpz(exp, e, j, NWORDS);
                extract_bignum_from_vec_to_mpz(mod, m, j, NWORDS);

                gmp_printf("init%d:\n\tbase = %Zx\n\texp = %Zx\n\tmod = %Zx\n",
                    j, base, exp, mod);
            }
        }

        // now we actually do the (vectorized) montgomery initialization
        // on our vector of random moduli.
        monty_init_vec(mtest, m, 0);

        if (verbose > 1)
        {
            for (j = 0; j < VECLEN; j++)
            {
                extract_bignum_from_vec_to_mpz(base, mtest->rhat, j, NWORDS);
                extract_bignum_from_vec_to_mpz(mod, mtest->one, j, NWORDS);

                gmp_printf("init%d:\n\trhat = %Zx\n\tone = %Zx\n\trho = %08x\n",
                    j, base, mod, mtest->vrho[j]);
            }
        }

        vecmulmod_ptr(b, mtest->rhat, b, m, s, mtest);      // monty rep
        vecmulmod_ptr(e, mtest->rhat, e, m, s, mtest);      // monty rep

        printf("thread %d starting %d iterations\n", tid, loc_iterations);

        // now do the calculation "b^e % m" a bunch of times
        for (i = 0; i < loc_iterations; i++)
        {
            if (verbose > 1)
            {
                for (j = 0; j < VECLEN; j++)
                {
                    extract_bignum_from_vec_to_mpz(base, b, j, NWORDS);

                    gmp_printf("monty(base%d) = %Zx\n", j, base);
                }
            }

            vecsqrmod_ptr(b, b, m, s, mtest);

            if (verbose > 1)
            {
                for (j = 0; j < VECLEN; j++)
                {
                    extract_bignum_from_vec_to_mpz(base, d, j, NWORDS);

                    gmp_printf("modexp%d = %Zx\n", j, base);
                }
            }

            // now verify each result
            if (do_verification)
            {
                vecmulmod_ptr(b, one, b, m, s, mtest);              // normal rep
                for (j = 0; j < VECLEN; j++)
                {
                    extract_bignum_from_vec_to_mpz(t1, d, j, NWORDS);
                    extract_bignum_from_vec_to_mpz(base, b, j, NWORDS);
                    extract_bignum_from_vec_to_mpz(exp, e, j, NWORDS);
                    extract_bignum_from_vec_to_mpz(mod, m, j, NWORDS);

                    mpz_powm(t2, base, exp, mod);

                    if (verbose)
                    {
                        gmp_printf("iteration %d lane %d:\n\tgmp = %Zx\n\ttest = %Zx\n",
                            i, j, t2, t1);
                    }

                    if (mpz_cmp(t1, t2) != 0)
                    {
                        gmp_printf("iteration %d error lane %d:\nbase = %Zx\nexp = %Zx\nmod = %Zx\ngmp = %Zx\ntest = %Zx\n",
                            i, j, base, exp, mod, t2, t1);
                        exit(1);
                    }

                }
            }
        }

        monty_free(mtest);

        if ((tid == 0) && (do_verification == 1))
            printf("verified %d x 16 vecModExp (all variable) results\n", loc_iterations);

        gettimeofday(&stopt, NULL);
        t_time = my_difftime(&startt, &stopt);
        elapsed_time[tid] = t_time;

        if (tid == 0)
            printf("Test with %d iterations took %1.4f seconds.\n", loc_iterations, t_time);

        mpz_clear(t1);
        mpz_clear(t2);
        mpz_clear(base);
        mpz_clear(mod);
        mpz_clear(exp);
        vecFree(m);
        vecFree(b);
        vecFree(d);
        vecFree(s);
        vecFree(e);
    }

    {
        int i;
        double sum = 0.0;
        double min_t = 9999999999.;
        double max_t = 0.;

        for (i = 0; i < threads; i++)
        {
            sum += elapsed_time[i];
            if (elapsed_time[i] < min_t)
                min_t = elapsed_time[i];
            if (elapsed_time[i] > max_t)
                max_t = elapsed_time[i];
        }

        printf("average elapsed time = %1.4f\n", sum / threads);
        printf("min elapsed time = %1.4f\n", min_t);
        printf("max elapsed time = %1.4f\n", max_t);
    }

    free(elapsed_time);
    free(LCG_STATE);

    printf("\n\n");

    return;
}

int main(int argc, char **argv)
{
    int threads;
    int do_verification = 1;
    int verbose = 0;

    if (argc < 2)
    {
        printf("usage: avx512_modexp $threads $do_verification\n");
        exit(1);
    }
    else if (argc == 3)
    {
        do_verification = atoi(argv[2]);
    }

    threads = atoi(argv[1]);

    printf("configured with MAXBITS = %d, DIGITBITS = %d, NUMWORDS = %d, VECLEN = %d\n",
        MAXBITS, DIGITBITS, NWORDS, VECLEN);
    printf("commencing modular exponentiation benchmarks\n"); fflush(stdout);

#ifdef BASE52
    vecmulmod_ptr = &vecmulmod52;
    vecsqrmod_ptr = &vecsqrmod52;
    montsetup_ptr = &vec_montgomery_setup52;
    vecmodexp_ptr = &vecmodexp52;
#else
    vecmulmod_ptr = &vecmulmod;
    vecsqrmod_ptr = &vecsqrmod;
    montsetup_ptr = &vec_montgomery_setup;
    vecmodexp_ptr = &vecmodexp;
#endif

    omp_set_num_threads(threads);    
    //vecmultest(do_verification, threads, verbose);
    //vecsqrtest(do_verification, threads, verbose);
    vecpmodtest(do_verification, threads, verbose);
    
    return 0;
}


