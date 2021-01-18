//	Copyright (c) 2021 by The Mayo Clinic, though its Special Purpose
//	 Processor Development Group (SPPDG). All Rights Reserved Worldwide.
//	 Licensed under the Apache License, Version 2.0 (the "License"); you may
//	 not use this file except in compliance with the License. You may obtain
//	 a copy of the License at http://www.apache.org/licenses/LICENSE-2.0.
//	Unless required by applicable law or agreed to in writing, software
//	 distributed under the License is distributed on an "AS IS" BASIS,
//	 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied,
//             including conditions of title, non-infringement, merchantability,
//            or fitness for a particular purpose
//	 See the License for the specific language governing permissions and
//	 limitations under the License.
//  This file is a snapshot of a work in progress, originated by Mayo
//   Clinic SPPDG. 

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

// this file defines a test harness for various modular exponentiation routines.
// We read command line options and execute the appropriate test(s).
#include "util.h"
#include "bigarith.h"
#include "pmod.h"
#include "monty_arith.h"
#include "x64_arith.h"
#include "gmp.h"

void mul_test(int iterations, int verbose, uint64_t lcg_state, pmod_t* pmod_state)
{
    int i, j, r;
    bignum* a, * b, * c, * n, * s;
    monty* mdata;
    uint64_t chksum = 0;

    struct timeval stopt;	// stop time of this job
    struct timeval startt;	// start time of this job
    double t_time = 0.;

    a = zInit();
    b = zInit();
    c = zInit();
    n = zInit();
    s = zInit();
    mdata = monty_alloc();

    for (j = 0; j < NWORDS; j++)
        a->data[j] = spRand64(&lcg_state);
    a->size = NWORDS;

    for (j = 0; j < NWORDS; j++)
        b->data[j] = spRand64(&lcg_state);
    b->size = NWORDS;

    for (j = 0; j < NWORDS; j++)
        n->data[j] = spRand64(&lcg_state);
    n->size = NWORDS;

    if ((n->data[0] & 1) == 0)
        n->data[0]++;

    // initialize the montgomery representation of this modulus.        
    monty_init(mdata, n, verbose);

    if (verbose > 0)
    {
        printf("initial a = "); zPrint(a); printf("\n");
        printf("initial b = "); zPrint(b); printf("\n");
        printf("initial n = "); zPrint(n); printf("\n");
    }

    to_monty(mdata, a);
    to_monty(mdata, b);

    gettimeofday(&startt, NULL);

    for (i = 0; i < iterations; i++)
    {
        if (verbose > 1)
        {
            printf("test %d:\n", i);
            printf("a = "); zPrint(a); printf("\n");
            printf("b = "); zPrint(b); printf("\n");
            printf("n = "); zPrint(n); printf("\n");
        }

        mul_ptr(mdata, a, b, a, n);
        chksum += a->data[0];

        if (verbose > 0)
        {
            printf("result: "); zPrint(a); printf("\n");
        }
    }

    gettimeofday(&stopt, NULL);
    t_time = my_difftime(&startt, &stopt);

    if (verbose > 0)
    {
        printf("final result: "); zPrint(a); printf("\n");
    }

    printf("final chksum: %lu\n", chksum);

    printf("%d mulredc tests took %.4f seconds\n", iterations, t_time);

    zFree(a);
    zFree(b);
    zFree(c);
    zFree(n);
    zFree(s);
    monty_free(mdata);

    return;
}

void sqr_test(int iterations, int verbose, uint64_t lcg_state, pmod_t* pmod_state)
{
    int i, j, r;
    bignum* a, * b, * c, * n, * s;
    monty* mdata;
    uint64_t chksum = 0;

    struct timeval stopt;	// stop time of this job
    struct timeval startt;	// start time of this job
    double t_time = 0.;

    a = zInit();
    b = zInit();
    c = zInit();
    n = zInit();
    s = zInit();
    mdata = monty_alloc();

    for (j = 0; j < NWORDS; j++)
        a->data[j] = spRand64(&lcg_state);
    a->size = NWORDS;

    for (j = 0; j < NWORDS; j++)
        n->data[j] = spRand64(&lcg_state);
    n->size = NWORDS;

    if ((n->data[0] & 1) == 0)
        n->data[0]++;

    // initialize the montgomery representation of this modulus.        
    monty_init(mdata, n, verbose);

    if (verbose > 0)
    {
        printf("initial a = "); zPrint(a); printf("\n");
        printf("initial n = "); zPrint(n); printf("\n");
    }

    to_monty(mdata, a);

    gettimeofday(&startt, NULL);

    for (i = 0; i < iterations; i++)
    {
        if (verbose > 1)
        {
            printf("test %d:\n", i);
            printf("a = "); zPrint(a); printf("\n");
            printf("n = "); zPrint(n); printf("\n");
        }

        sqr_ptr(mdata, a, a, n);
        chksum += a->data[0];

        if (verbose > 1)
        {
            printf("result: "); zPrint(a); printf("\n");
        }
    }

    gettimeofday(&stopt, NULL);
    t_time = my_difftime(&startt, &stopt);

    if (verbose > 0)
    {
        printf("final result: "); zPrint(a); printf("\n");
    }

    printf("final chksum: %lu\n", chksum);
    printf("%d sqrredc tests took %.4f seconds\n", iterations, t_time);

    zFree(a);
    zFree(b);
    zFree(c);
    zFree(n);
    zFree(s);
    monty_free(mdata);

    return;
}

void monty_test(int iterations, int verbose, uint64_t lcg_state, pmod_t *pmod_state)
{
    int i, j, r;
    bignum *a, *b, *c, *n, *s;
    monty *mdata;

    struct timeval stopt;	// stop time of this job
    struct timeval startt;	// start time of this job
    double t_time = 0.;

    a = zInit();
    b = zInit();
    c = zInit();
    n = zInit();
    s = zInit();
    mdata = monty_alloc();    

    gettimeofday(&startt, NULL);

    for (i = 0; i < iterations; i++)
    {
        for (j = 0; j < NWORDS; j++)
            a->data[j] = spRand64(&lcg_state);
        a->size = NWORDS;

        for (j = 0; j < NWORDS; j++)
            b->data[j] = spRand64(&lcg_state);
        b->size = NWORDS;

        for (j = 0; j < NWORDS; j++)
            n->data[j] = spRand64(&lcg_state);
        n->size = NWORDS;

        if ((n->data[0] & 1) == 0)
            n->data[0]++;        

        // initialize the montgomery representation of this modulus.        
        monty_init(mdata, n, verbose);

        if (verbose)
        {
            printf("test %d:\n", i);
            printf("a = "); zPrint(a); printf("\n");
            printf("b = "); zPrint(b); printf("\n");
            printf("n = "); zPrint(n); printf("\n");
        }

        to_monty(mdata, a);

        lroddwin_powm(pmod_state, mdata, c, a, b, n, s);

        if (verbose)
        {
            printf("result: "); zPrint(c); printf("\n");
        }
    }

    gettimeofday(&stopt, NULL);
    t_time = my_difftime(&startt, &stopt);

    printf("%d powm tests took %.4f seconds\n", iterations, t_time);

    zFree(a);
    zFree(b);
    zFree(c);
    zFree(n);
    zFree(s);
    monty_free(mdata);

    return;
}

int main(int argc, char **argv)
{
    struct timeval stopt;	// stop time of this job
    struct timeval startt;	// start time of this job
    double t_time = 0.;
    int iterations = 1000;
    int seed;
    uint64_t *lcg_state;
    pmod_t *pmod_state;
    int verbose = 0;
    int taskid, numtasks;

    if (argc > 1)
    {
        iterations = atoi(argv[1]);
    }

    if (argc > 2)
    {
        verbose = atoi(argv[2]);
    }

    if (argc > 3)
    {
        seed = atoi(argv[3]);
    }
    else
    {
        gettimeofday(&startt, NULL);
        seed = hash64((startt.tv_usec));
    }


    lcg_state = (uint64_t *)malloc(1 * sizeof(uint64_t));
    lcg_state[0] = hash64((seed));
    pmod_state = (pmod_t *)malloc(sizeof(pmod_t));
    pmodlib_init(pmod_state);

    printf("commencing benchmarks with MAXBITS = %d, NWORDS = %d\n",
        MAXBITS, NWORDS);

    // configure benchmark tests
    int do_pmod_tests = 1;
    int do_mulsqr_tests = 1;

    int bench_sos = 1;
    int bench_fios = 1;
    int bench_fips = 1;
    int bench_cios = 1;
    int bench_bps = 1;
    int bench_gmp = 1;

    if (do_mulsqr_tests)
    {
        if (bench_sos)
        {
            printf("commencing %d mulredc iterations using mulmod_sos\n", iterations);
            mul_ptr = &mulmod_sos;
            mul_test(iterations, verbose, lcg_state[0], pmod_state);
        }

        if (bench_cios)
        {
            printf("commencing %d mulredc iterations using mulmod_cios\n", iterations);
            mul_ptr = &mulmod_cios;
            mul_test(iterations, verbose, lcg_state[0], pmod_state);
        }

        if (bench_bps)
        {
            printf("commencing %d mulredc iterations using mulmod_bps\n", iterations);
            mul_ptr = &mulmod_bps;
            mul_test(iterations, verbose, lcg_state[0], pmod_state);
        }

        if (bench_fios)
        {
            printf("commencing %d mulredc iterations using mulmod_fios\n", iterations);
            mul_ptr = &mulmod_fios;
            mul_test(iterations, verbose, lcg_state[0], pmod_state);
        }

        if (bench_fips)
        {
            printf("commencing %d mulredc iterations using mulmod_fips\n", iterations);
            mul_ptr = &mulmod_fips;
            mul_test(iterations, verbose, lcg_state[0], pmod_state);
        }

        if (bench_sos)
        {
            printf("commencing %d sqrredc iterations using sqrmod_sos\n", iterations);
            sqr_ptr = &sqrmod_sos;
            sqr_test(iterations, verbose, lcg_state[0], pmod_state);
        }

        if (bench_cios)
        {
            printf("commencing %d sqrredc iterations using sqrmod_cios\n", iterations);
            sqr_ptr = &sqrmod_cios;
            sqr_test(iterations, verbose, lcg_state[0], pmod_state);
        }

        if (bench_bps)
        {
            printf("commencing %d sqrredc iterations using sqrmod_bps\n", iterations);
            sqr_ptr = &sqrmod_bps;
            sqr_test(iterations, verbose, lcg_state[0], pmod_state);
        }

        if (bench_fios)
        {
            printf("commencing %d sqrredc iterations using sqrmod_fios\n", iterations);
            sqr_ptr = &sqrmod_fios;
            sqr_test(iterations, verbose, lcg_state[0], pmod_state);
        }

        if (bench_fips)
        {
            printf("commencing %d sqrredc iterations using sqrmod_fips\n", iterations);
            sqr_ptr = &sqrmod_fips;
            sqr_test(iterations, verbose, lcg_state[0], pmod_state);
        }

        // gmp SOS mul
        if (bench_gmp)
        {
            mpz_t a, b, n, t, nhat, r, u;
            int i, j, k;
            struct timeval stopt;	// stop time of this job
            struct timeval startt;	// start time of this job
            double t_time = 0.;
            int numbits = MAXBITS;
            uint64_t lcg = lcg_state[0];
            uint64_t chksum = 0;

            mpz_init(a);
            mpz_init(b);
            mpz_init(nhat);
            mpz_init(r);
            mpz_init(n);
            mpz_init(t);
            mpz_init(u);

            gettimeofday(&startt, NULL);

            mpz_set_ui(a, 0);
            for (j = 0; j < NWORDS; j++)
            {
                uint64_t x = spRand64(&lcg);
                mpz_set_ui(t, x);
                mpz_mul_2exp(t, t, 64 * j);
                mpz_add(a, a, t);
            }

            mpz_set_ui(b, 0);
            for (j = 0; j < NWORDS; j++)
            {
                uint64_t x = spRand64(&lcg);
                mpz_set_ui(t, x);
                mpz_mul_2exp(t, t, 64 * j);
                mpz_add(b, b, t);
            }

            mpz_set_ui(n, 0);
            for (j = 0; j < NWORDS; j++)
            {
                uint64_t x = spRand64(&lcg);
                mpz_set_ui(t, x);
                mpz_mul_2exp(t, t, 64 * j);
                mpz_add(n, n, t);
            }

            if ((mpz_get_ui(n) & 1) == 0)
            {
                mpz_add_ui(n, n, 1);
            }

            printf("commencing %d mulredc iterations using gmp \n", iterations);
            
            if (verbose > 0)
            {
                gmp_printf("initial a: %Zx\n", a);
                gmp_printf("initial b: %Zx\n", b);
                gmp_printf("initial n: %Zx\n", n);
            }

            // monty setup
            mpz_set_ui(r, 1);
            mpz_mul_2exp(r, r, MAXBITS);
            mpz_invert(nhat, n, r);
            mpz_sub(nhat, r, nhat);
            mpz_mul(a, r, a);
            mpz_tdiv_r(a, a, n);
            mpz_mul(b, r, b);
            mpz_tdiv_r(b, b, n);

            gettimeofday(&startt, NULL);
            for (k = 0; k < iterations; k++)
            {
                mpz_mul(t, a, b);
                mpz_tdiv_r_2exp(a, t, MAXBITS);
                mpz_mul(u, a, nhat);
                mpz_tdiv_r_2exp(u, u, MAXBITS);
                mpz_mul(a, u, n);
                mpz_add(a, t, a);
                mpz_tdiv_q_2exp(a, a, MAXBITS);
                if (mpz_sizeinbase(a, 2) > MAXBITS)
                    mpz_sub(a, a, n);

                if (verbose > 0)
                {
                    gmp_printf("result: %Zx\n", a);
                }
                chksum += a->_mp_d[0];
            }

            gettimeofday(&stopt, NULL);
            t_time = my_difftime(&startt, &stopt);

            if (verbose > 0)
            {
                gmp_printf("final result: %Zx\n", a);
            }

            printf("final chksum: %lu\n", chksum);
            printf("%d mulredc tests took %.4f seconds\n", iterations, t_time);

            mpz_clear(a);
            mpz_clear(b);
            mpz_clear(r);
            mpz_clear(nhat);
            mpz_clear(n);
            mpz_clear(t);
            mpz_clear(u);
        }

        // gmp SOS sqr
        if (bench_gmp)
        {
            mpz_t a, b, n, t, nhat, r, u;
            int i, j, k;
            struct timeval stopt;	// stop time of this job
            struct timeval startt;	// start time of this job
            double t_time = 0.;
            int numbits = MAXBITS;
            uint64_t lcg = lcg_state[0];
            uint64_t chksum = 0;

            mpz_init(a);
            mpz_init(b);
            mpz_init(nhat);
            mpz_init(r);
            mpz_init(n);
            mpz_init(t);
            mpz_init(u);

            gettimeofday(&startt, NULL);

            mpz_set_ui(a, 0);
            for (j = 0; j < NWORDS; j++)
            {
                uint64_t x = spRand64(&lcg);
                mpz_set_ui(t, x);
                mpz_mul_2exp(t, t, 64 * j);
                mpz_add(a, a, t);
            }

            mpz_set_ui(n, 0);
            for (j = 0; j < NWORDS; j++)
            {
                uint64_t x = spRand64(&lcg);
                mpz_set_ui(t, x);
                mpz_mul_2exp(t, t, 64 * j);
                mpz_add(n, n, t);
            }

            if ((mpz_get_ui(n) & 1) == 0)
            {
                mpz_add_ui(n, n, 1);
            }

            printf("commencing %d sqrredc iterations using gmp \n", iterations);

            if (verbose > 0)
            {
                gmp_printf("initial a: %Zx\n", a);
                gmp_printf("initial n: %Zx\n", n);
            }

            // monty setup
            mpz_set_ui(r, 1);
            mpz_mul_2exp(r, r, MAXBITS);
            mpz_invert(nhat, n, r);
            mpz_sub(nhat, r, nhat);
            mpz_mul(a, r, a);
            mpz_tdiv_r(a, a, n);

            gettimeofday(&startt, NULL);
            for (k = 0; k < iterations; k++)
            {
                mpz_mul(t, a, a);
                mpz_tdiv_r_2exp(a, t, MAXBITS);
                mpz_mul(u, a, nhat);
                mpz_tdiv_r_2exp(u, u, MAXBITS);
                mpz_mul(a, u, n);
                mpz_add(a, t, a);
                mpz_tdiv_q_2exp(a, a, MAXBITS);
                if (mpz_sizeinbase(a, 2) > MAXBITS)
                    mpz_sub(a, a, n);

                chksum += mpz_get_ui(a);
                if (verbose > 0)
                {
                    gmp_printf("result: %Zx\n", a);
                }
            }

            gettimeofday(&stopt, NULL);
            t_time = my_difftime(&startt, &stopt);

            if (verbose > 0)
            {
                gmp_printf("final result: %Zx\n", a);
            }

            printf("final chksum: %lu\n", chksum);
            printf("%d sqrredc tests took %.4f seconds\n", iterations, t_time);

            mpz_clear(a);
            mpz_clear(b);
            mpz_clear(r);
            mpz_clear(nhat);
            mpz_clear(n);
            mpz_clear(t);
            mpz_clear(u);
        }
    }

    if (do_mulsqr_tests)
        iterations /= 10000;

    if (do_pmod_tests)
    {
        if (bench_sos)
        {
            printf("commencing %d powm iterations using mulmod_sos\n", iterations);
            mul_ptr = &mulmod_sos;
            sqr_ptr = &sqrmod_sos_mul;
            monty_test(iterations, verbose, lcg_state[0], pmod_state);
        }

        if (bench_sos)
        {
            printf("commencing %d powm iterations using mulmod_sos and sqrmod_sos\n", iterations);
            mul_ptr = &mulmod_sos;
            sqr_ptr = &sqrmod_sos;
            monty_test(iterations, verbose, lcg_state[0], pmod_state);
        }

        if (bench_cios)
        {
            printf("commencing %d powm iterations using mulmod_cios\n", iterations);
            mul_ptr = &mulmod_cios;
            sqr_ptr = &sqrmod_cios_mul;
            monty_test(iterations, verbose, lcg_state[0], pmod_state);
        }

        if (bench_cios)
        {
            printf("commencing %d powm iterations using mulmod_cios and sqrmod_cios\n", iterations);
            mul_ptr = &mulmod_cios;
            sqr_ptr = &sqrmod_cios;
            monty_test(iterations, verbose, lcg_state[0], pmod_state);
        }

        if (bench_bps)
        {
            printf("commencing %d powm iterations using mulmod_bps\n", iterations);
            mul_ptr = &mulmod_bps;
            sqr_ptr = &sqrmod_bps_mul;
            monty_test(iterations, verbose, lcg_state[0], pmod_state);
        }

        if (bench_bps)
        {
            printf("commencing %d powm iterations using mulmod_bps and sqrmod_bps\n", iterations);
            mul_ptr = &mulmod_bps;
            sqr_ptr = &sqrmod_bps;
            monty_test(iterations, verbose, lcg_state[0], pmod_state);
        }

        if (bench_fios)
        {
            printf("commencing %d powm iterations using mulmod_fios\n", iterations);
            mul_ptr = &mulmod_fios;
            sqr_ptr = &sqrmod_fios_mul;
            monty_test(iterations, verbose, lcg_state[0], pmod_state);
        }

        if (bench_fios)
        {
            printf("commencing %d powm iterations using mulmod_fios and sqrmod_fios\n", iterations);
            mul_ptr = &mulmod_fios;
            sqr_ptr = &sqrmod_fios;
            monty_test(iterations, verbose, lcg_state[0], pmod_state);

        }

        if (bench_fips)
        {
            printf("commencing %d powm iterations using mulmod_fips\n", iterations);
            mul_ptr = &mulmod_fips;
            sqr_ptr = &sqrmod_fips_mul;
            monty_test(iterations, verbose, lcg_state[0], pmod_state);
        }

        if (bench_fips)
        {
            printf("commencing %d powm iterations using mulmod_fips and sqrmod_fips\n", iterations);
            mul_ptr = &mulmod_fips;
            sqr_ptr = &sqrmod_fips;
            monty_test(iterations, verbose, lcg_state[0], pmod_state);
        }

        // gmp comparison.  These results won't match the ones above because
        // we use a different RNG (the builtin gmp RNG).
        if (bench_gmp)
        {
            mpz_t a, b, n, t, aa, bb;
            int i, j, k;
            struct timeval stopt;	// stop time of this job
            struct timeval startt;	// start time of this job
            double t_time = 0.;
            int numbits = MAXBITS;
            gmp_randstate_t gmp_randstate;

            mpz_init(a);
            mpz_init(b);
            mpz_init(aa);
            mpz_init(bb);
            mpz_init(n);
            mpz_init(t);

            gettimeofday(&startt, NULL);
            srand(42); // lcg_state[0]);
            gmp_randinit_default(gmp_randstate);
            gmp_randseed_ui(gmp_randstate, rand());

            printf("commencing %d powm iterations using gmp powm\n", iterations);
            mpz_urandomb(n, gmp_randstate, numbits);
            if (mpz_even_p(n))
                mpz_add_ui(n, n, 1);

            gettimeofday(&startt, NULL);
            for (k = 0; k < iterations; k++)
            {
                mpz_urandomb(a, gmp_randstate, numbits);
                mpz_urandomb(b, gmp_randstate, numbits);

                mpz_tdiv_r(a, a, n);
                mpz_tdiv_r(b, b, n);

                if (verbose)
                {
                    printf("test %d:\n", k);
                    printf("a = "); gmp_printf("%Zx\n", a);
                    printf("b = "); gmp_printf("%Zx\n", b);
                    printf("n = "); gmp_printf("%Zx\n", n);
                }

                mpz_set(aa, a);
                mpz_set(bb, b);

                mpz_powm(a, aa, bb, n);

                if (verbose)
                {
                    gmp_printf("result: %Zx\n", a);
                }
            }

            gettimeofday(&stopt, NULL);
            t_time = my_difftime(&startt, &stopt);

            printf("%d powm tests took %.4f seconds\n", iterations, t_time);

            mpz_clear(a);
            mpz_clear(b);
            mpz_clear(aa);
            mpz_clear(bb);
            mpz_clear(n);
            mpz_clear(t);
        }
    }

    free(lcg_state);
    pmodlib_free(pmod_state);
    free(pmod_state);
  
    return 0;
}
