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
#include "pmod.h"


void pmodlib_init(pmod_t *pmod_state)
{
    int i;

    // accomodate window sizes up to 8
    pmod_state->libpmod_gwin = (bignum **)malloc((1 << MAX_WINSIZE) * sizeof(bignum *));
    pmod_state->libpmod_gwin[0] = zInit();

    for (i = 1; i < (1 << MAX_WINSIZE); i++)
    {
        pmod_state->libpmod_gwin[i] = zInit();
    }

    return;
}

void pmodlib_free(pmod_t *pmod_state)
{
    int i;

    for (i = 0; i < (1 << MAX_WINSIZE); i++)
    {
        zFree(pmod_state->libpmod_gwin[i]);
    }
    free(pmod_state->libpmod_gwin);

    return;
}

void lr_powm(pmod_t *pmod_state, monty *mdata, bignum *c, bignum *a, bignum *b, bignum *n, bignum *s)
{
    bignum *acc;
    int i;
    int j;

    acc = zInit();
    zCopy(mdata->one, acc);

    for (i = NWORDS - 1; i >= 0; i--)
    {
        for (j = 63; j >= 0; j--)
        {
            sqr_ptr(mdata, acc, acc, s);
            if (b->data[i] & (1ULL << j))
                mul_ptr(mdata, acc, a, acc, s);
        }
    }

    memset(mdata->mtmp1->data, 0, 2 * NWORDS * sizeof(base_t));
    zSet1(mdata->mtmp1, 1);
    mul_ptr(mdata, acc, mdata->mtmp1, c, s);

    // final check to ensure c < N
    i = 1;
    for (j = NWORDS - 1; j >= 0; j--)
    {
        if (c->data[j] > mdata->n->data[j])
            break;

        if (c->data[j] < mdata->n->data[j])
        {
            i = 0;
            break;
        }
    }

    if (i)
    {
        mpSub(c->data, mdata->n->data, c->data, NWORDS);
    }

    zFree(acc);
    c->size = NWORDS;

    return;
}

int get_winsize(void)
{
    // the window size is based on minimizing the total number of multiplications
    // in the windowed exponentiation.  experiments show that this is best;
    // the growing size of the table doesn't change the calculus, at least
    // on the KNL.
    int size;
    int muls;
    int minmuls = 99999999;
    int minsize = 4;

    for (size = 2; size <= 8; size++)
    {
        muls = (MAXBITS / size) + (1 << size);
        if (muls < minmuls)
        {
            minmuls = muls;
            minsize = size;
        }
    }

    return minsize;
}

int get_bitwin(bignum *b, int bitloc, int winsize, int winmask)
{
    int bstr;
    int bitstart = (bitloc - winsize + 1);
    int word = bitloc / 64;
    int word2 = bitstart / 64;

    bitstart = bitstart % 64;

    if (word == word2)
    {
        bstr = (b->data[word] >> bitstart) & winmask;
    }
    else
    {
        int upperbits = (bitloc % 64) + 1;

        bstr = (b->data[word2] >> bitstart);
        bstr |= ((b->data[word]) << (winsize - upperbits));
        bstr &= winmask;
    }

    return bstr;
}

int get_oddbitwin(bignum *b, int bitloc, int winsize, int winmask, int *m)
{
    int bstr;
    int bitstart = (bitloc - winsize + 1);
    int word = bitloc / 64;
    int word2 = bitstart / 64;

    bitstart = bitstart % 64;

    if (word == word2)
    {
        bstr = (b->data[word] >> bitstart) & winmask;
    }
    else
    {
        int upperbits = (bitloc % 64) + 1;

        bstr = (b->data[word2] >> bitstart);
        bstr |= ((b->data[word]) << (winsize - upperbits));
        bstr &= winmask;
    }

    *m = 0;
    while ((bstr & 1) == 0)
    {
        if (bstr == 0)
            break;

        (*m)++;
        bstr >>= 1;
    }

    return bstr;
}

void lrwin_powm(pmod_t *pmod_state, monty *mdata, bignum *c, bignum *a, bignum *b, bignum *n, bignum *s)
{
    bignum *acc;
    int i, j, bit = MAXBITS - 1;
    int k = get_winsize();
    bignum **g = pmod_state->libpmod_gwin;             // storage for windowed method precomputation
    int mask;
    int bstr;

    mask = 0;
    for (j = 0; j < k; j++)
    {
        mask = (mask << 1) | 1;
    }

    acc = zInit();
    zCopy(mdata->one, acc);

    // precomputations, b^i for 0 <= i < 2^k
    memcpy(g[1]->data, a->data, NWORDS * sizeof(base_t));
    for (i = 2; i < (1 << k); i++)
    {
        mul_ptr(mdata, g[i - 1], a, g[i], s);
    }

    // L-R windowed exponentiation.  Scan the exponent bit-vector
    // backward instead of flipping and shifting it.
    while (bit >= 0)
    {
        if (bit < k)
        {
            // grab the last bits of the exponent.
            // accommodates exponent lengths not divisible
            // by the window size
            mask = 0x0;
            for (j = 0; j < (bit + 1); j++)
            {
                sqr_ptr(mdata, acc, acc, s);
                mask = (mask << 1) | 1;
            }

            bstr = b->data[0] & mask;
        }
        else
        {
            // grab the next k bits of the exponent.
            bstr = get_bitwin(b, bit, k, mask);
            for (j = 0; j < k; j++)
            {
                sqr_ptr(mdata, acc, acc, s);
            }
        }

        if (bstr > 0)
            mul_ptr(mdata, acc, g[bstr], acc, s);

        bit -= k;

    }

    memset(mdata->mtmp1->data, 0, 2 * NWORDS * sizeof(base_t));
    zSet1(mdata->mtmp1, 1);
    mul_ptr(mdata, acc, mdata->mtmp1, c, s);

    // final check to ensure c < N
    i = 1;
    for (j = NWORDS - 1; j >= 0; j--)
    {
        if (c->data[j] > mdata->n->data[j])
            break;

        if (c->data[j] < mdata->n->data[j])
        {
            i = 0;
            break;
        }
    }

    if (i)
    {
        mpSub(c->data, mdata->n->data, c->data, NWORDS);
    }
    c->size = NWORDS;

    zFree(acc);

    return;
}

void lroddwin_powm(pmod_t *pmod_state, monty *mdata, bignum *c, bignum *a, bignum *b, bignum *n, bignum *s)
{
    bignum *acc;
    int i, j, bit = MAXBITS - 1;
    int k = get_winsize();
    bignum **g = pmod_state->libpmod_gwin;             // storage for windowed method precomputation
    int mask;
    int bstr;
    int m;

    mask = 0;
    for (j = 0; j < k; j++)
    {
        mask = (mask << 1) | 1;
    }

    acc = zInit();
    zCopy(mdata->one, acc);

    // precomputations, b^i for 0 <= i < 2^k, i odd (except i = 2).
    // half the setup cost for minimal extra overhead while scanning
    // the exponent vector... not as secure because order of
    // operations (squaring/multiply) depends on exponent bits.
    memcpy(g[1]->data, a->data, NWORDS * sizeof(base_t));
    mul_ptr(mdata, g[1], a, g[2], s);
    //printf("g[%d] ", 2); zPrint(g[2]); printf("\n");
    for (i = 3; i < (1 << k); i += 2)
    {
        mul_ptr(mdata, g[i - 2], g[2], g[i], s);
        //printf("g[%d] ", i); zPrint(g[i]); printf("\n");
    }

    //printf("acc init "); zPrint(acc); printf("\n");

    // L-R windowed exponentiation.  Scan the exponent bit-vector
    // backward instead of flipping and shifting it.
    while (bit >= 0)
    {
        if (bit < (k- 1))
        {
            // grab the last bits of the exponent.
            // accommodates exponent lengths not divisible
            // by the window size
            mask = 0x0;
            k = (bit + 1);
            for (j = 0; j < k; j++)
            {
                mask = (mask << 1) | 1;
            }
        }

        // grab the next k bits of the exponent.
        bstr = get_oddbitwin(b, bit, k, mask, &m);
        for (j = 0; j < (k - m); j++)
        {
            sqr_ptr(mdata, acc, acc, s);
            //printf("sqr bit %03d ", bit); zPrint(acc); printf("\n");
        }

        if (bstr > 0)
        {
            mul_ptr(mdata, acc, g[bstr], acc, s);
            //printf("mul bit %03d ", bit); zPrint(acc); printf("\n");
        }

        for (j = 0; j < m; j++)
        {
            sqr_ptr(mdata, acc, acc, s);
            //printf("sqr bit %03d ", bit); zPrint(acc); printf("\n");
        }

        bit -= k;
    }

    memset(mdata->mtmp1->data, 0, 2*NWORDS * sizeof(base_t));
    zSet1(mdata->mtmp1, 1);
    mul_ptr(mdata, acc, mdata->mtmp1, c, s);

    // final check to ensure c < N
    i = 1;
    for (j = NWORDS - 1; j >= 0; j--)
    {
        if (c->data[j] > mdata->n->data[j])
            break;

        if (c->data[j] < mdata->n->data[j])
        {
            i = 0;
            break;
        }
    }

    if (i)
    {
        mpSub(c->data, mdata->n->data, c->data, NWORDS);
    }
    c->size = NWORDS;

    zFree(acc);

    return;
}

