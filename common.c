#include "vecarith.h"

void(*vecmulmod_ptr)(bignum*, bignum*, bignum*, bignum*, bignum*, monty*);
void(*vecsqrmod_ptr)(bignum*, bignum*, bignum*, bignum*, monty*);
int(*montsetup_ptr)(bignum*, bignum*, bignum*, base_t*);
void(*vecmodexp_ptr)(bignum*, bignum*, bignum*, bignum*, bignum*, bignum*, monty* m);

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
        muls = (NWORDS * DIGITBITS / size) + (1 << size);
        if (muls < minmuls)
        {
            minmuls = muls;
            minsize = size;
        }
    }

    return minsize;
}

int get_bitwin(bignum *e, int bitloc, int winsize, int lane, int winmask)
{
    int bstr;
    int bitstart = (bitloc - winsize + 1);
    int word = bitloc / DIGITBITS;
    int word2 = bitstart / DIGITBITS;

    bitstart = bitstart % DIGITBITS;

    if (word == word2)
    {
        bstr = (e->data[lane + word * VECLEN] >> bitstart) & winmask;
    }
    else
    {
        int upperbits = (bitloc % DIGITBITS) + 1;

        bstr = (e->data[lane + word2 * VECLEN] >> bitstart);
        bstr |= ((e->data[lane + word * VECLEN]) << (winsize - upperbits));
        bstr &= winmask;
    }

    return bstr;
}


bignum * vecInit(void)
{
    int i;
    size_t sz = VECLEN * (2 * NWORDS + 4);
    bignum *n;
    n = (bignum *)malloc(sizeof(bignum));

    n->data = (base_t *)xmalloc_align(sz * sizeof(base_t));
    if (n->data == NULL)
    {
        printf("could not allocate memory\n");
        exit(2);
    }

    for (i = 0; i < sz; i++)
    {
        n->data[i] = 0;
    }
    n->size = 1;

    return n;
}

void vecCopy(bignum * src, bignum * dest)
{
    //physically copy the digits of u into the digits of v
    int su = VECLEN * (2 * NWORDS + 1);

    memcpy(dest->data, src->data, su * sizeof(base_t));
    dest->size = src->size; // = NWORDS;
    return;
}

void vecCopyn(bignum * src, bignum * dest, int size)
{
    //physically copy the digits of u into the digits of v
    int su = VECLEN * size;

    memcpy(dest->data, src->data, su * sizeof(base_t));
    dest->size = size;
    return;
}

void vecClear(bignum *n)
{
    memset(n->data, 0, VECLEN*(2 * NWORDS + 1) * sizeof(base_t));
    return;
}

void vecFree(bignum *n)
{
    align_free(n->data);
    free(n);
}

void copy_vec_lane(bignum *src, bignum *dest, int num, int size)
{
    int j;

    for (j = 0; j < size; j++)
    {
        dest->data[num + j * VECLEN] = src->data[num + j * VECLEN];
    }

    return;
}

void monty_init_vec(monty *mdata, bignum * n, int verbose)
{
    int j;
    // for a input modulus n, initialize constants for 
    // montogomery representation
    // this assumes that n is relatively prime to 2, i.e. is odd.	
    // In this version we assume the input monty structure has
    // already been allocated and we just perform the calculations.

    if (verbose)
        printf("initializing montgomery representation\n");

    memset(mdata->n->data, 0, (2 * NWORDS * VECLEN + 1) * sizeof(base_t));
    memset(mdata->r->data, 0, (2 * NWORDS * VECLEN + 1) * sizeof(base_t));
    memset(mdata->rhat->data, 0, (2 * NWORDS * VECLEN + 1) * sizeof(base_t));
    memset(mdata->one->data, 0, (2 * NWORDS * VECLEN + 1) * sizeof(base_t));
    memset(mdata->vrho, 0, VECLEN * sizeof(base_t));

    vecCopy(n, mdata->n);
    montsetup_ptr(mdata->n, mdata->r, mdata->rhat, mdata->vrho);

    for (j = 0; j < VECLEN; j++)
    {
        mdata->one->data[j] = 1;
    }

    vecmulmod_ptr(mdata->one, mdata->rhat, mdata->one, n, mdata->mtmp1, mdata);    // monty rep
    vecCopyn(mdata->one, mdata->g[0], NWORDS);

    return;

}

monty* monty_alloc(void)
{
    int i;
    monty *mdata = (monty *)malloc(sizeof(monty));

    mdata->r = vecInit();
    mdata->n = vecInit();
    mdata->nhat = vecInit();
    mdata->vnhat = vecInit();
    mdata->rhat = vecInit();
    mdata->rmask = vecInit();
    mdata->one = vecInit();
    mdata->mtmp1 = vecInit();
    mdata->mtmp2 = vecInit();
    mdata->mtmp3 = vecInit();
    mdata->mtmp4 = vecInit();

    mdata->g = (bignum **)malloc((1 << MAX_WINSIZE) * sizeof(bignum *));
    mdata->g[0] = vecInit();

    for (i = 1; i < (1 << MAX_WINSIZE); i++)
    {
        mdata->g[i] = vecInit();
    }

    mdata->vrho = (base_t *)xmalloc_align(VECLEN * sizeof(base_t));

    return mdata;
}

void monty_free(monty *mdata)
{
    int i;

    vecFree(mdata->mtmp1);
    vecFree(mdata->mtmp2);
    vecFree(mdata->mtmp3);
    vecFree(mdata->mtmp4);
    vecFree(mdata->rhat);
    vecFree(mdata->one);
    vecFree(mdata->n);
    vecFree(mdata->nhat);
    vecFree(mdata->r);
    align_free(mdata->vrho);

    for (i = 0; i < (1 << MAX_WINSIZE); i++)
    {
        vecFree(mdata->g[i]);
    }
    free(mdata->g);

    return;
}

