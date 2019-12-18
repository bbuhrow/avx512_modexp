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

    mpz_set_ui(dest, 0);
    for (j = sz - 1; j >= 0; j--)
    {
        mpz_mul_2exp(dest, dest, sizeof(base_t) * 8);
        mpz_add_ui(dest, dest, vec_src->data[num + j * VECLEN]);
    }

    return;
}

void vecpmodtest(int do_verification)
{
    // test the pmod by comparing all results to those computed using
    // validated scalar code.
    int i, j;
    int iterations = 100000;
    int threads = 1;
    int verbose = 0;

    printf("commencing test: all variable (random)\n");
#pragma omp parallel num_threads(threads)
    {
        // timing variables
        struct timeval stopt;	// stop time of this job
        struct timeval startt;	// start time of this job
        double t_time;

        mpz_t base, exp, mod, t1, t2;
        gmp_randstate_t rng_state;
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
        gmp_randinit_default(rng_state);
        gmp_randseed_ui(rng_state, tid);

        // attempt to scale the number of iterations with input size
        // so this doesn't take forever.
        loc_iterations = iterations * 2 / (NWORDS * DIGITBITS);

        if (MAXBITS >= 4096)
            loc_iterations *= 1; 
        else if (MAXBITS >= 2048)
            loc_iterations *= 10;
        else if (MAXBITS >= 1024)
            loc_iterations *= 5;
        else if (MAXBITS >= 512)
            loc_iterations *= 10;
        else if (MAXBITS >= 256)
            loc_iterations *= 25;
        else
            loc_iterations *= 100;

        mtest = monty_alloc();
        gettimeofday(&startt, NULL);
        
        for (j = 0; j < VECLEN; j++)
        {
            one->data[j] = 1;
        }

        // now do the calculation "b^e % m" a bunch of times
        for (i = 0; i < loc_iterations; i++)
        {
            mpz_urandomb(mod, rng_state, MAXBITS * VECLEN);
            mpz_urandomb(exp, rng_state, MAXBITS * VECLEN);
            mpz_urandomb(base, rng_state, MAXBITS * VECLEN);

            memset(m->data, 0, MAXBITS * 2 * VECLEN / 8);
            memset(e->data, 0, MAXBITS * 2 * VECLEN / 8);
            memset(b->data, 0, MAXBITS * 2 * VECLEN / 8);
            memcpy(m->data, mod->_mp_d, MAXBITS * VECLEN / 8);
            memcpy(e->data, exp->_mp_d, MAXBITS * VECLEN / 8);
            memcpy(b->data, base->_mp_d, MAXBITS * VECLEN / 8);
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

            vecmulmod(b, mtest->rhat, b, m, s, mtest);      // monty rep

            if (verbose > 1)
            {
                for (j = 0; j < VECLEN; j++)
                {
                    extract_bignum_from_vec_to_mpz(base, b, j, NWORDS);

                    gmp_printf("monty(base%d) = %Zx\n", j, base);
                }
            }

            vecmodexp(d, b, e, m, s, mtest->one, mtest);    // powm
            vecmulmod(d, one, d, m, s, mtest);              // normal rep

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
                vecmulmod(b, one, b, m, s, mtest);              // normal rep
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

    printf("\n\n");

    return;
}


int main(int argc, char **argv)
{
    int i;
    int threads;
    int do_verification = 1;

    if (argc < 2)
    {
        printf("usage: phi_modexp $threads $do_verification\n");
        exit(1);
    }
    else if (argc == 3)
    {
        do_verification = atoi(argv[2]);
    }

    threads = 1;
    if (argc == 2)
        threads = atoi(argv[1]);

    printf("commencing modular exponentiation benchmarks using %d threads and window size %d "
        "for MAXBITS = %d\n", threads, get_winsize(), MAXBITS); fflush(stdout);

    omp_set_num_threads(threads);    
    vecpmodtest(do_verification);

    return 0;
}


