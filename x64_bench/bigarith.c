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

#include "bigarith.h"

void zSet1(bignum *n, base_t d)
{
    n->data[0] = d;
    n->size = 1;
    return;
}

int zBits(bignum * n)
{
	if (n->size == 1)
		return spBits(n->data[0]);
	else
		return DIGITBITS*(n->size-1) + spBits(n->data[n->size-1]);
}

base_t spBits(base_t n)
{
	int i = 0;
	while (n != 0)
	{
		n >>= 1;
		i++;
	}
	return i;
}

int ndigits_1(base_t n)
{
	int i=0;
	while (n != 0)
	{
		n /= 10;
		i++;
	}
	if (i==0)
		i++;
	return i;
}

base_t spGCD(base_t x, base_t y)
{
	base_t a,b,c;
	a=x; b=y;
	while (b != 0)
	{
		c=a%b;
		a=b;
		b=c;
	}
	return a;
}

void sp2big(base_t src, bignum * dest)
{
	dest->data[0] = src;
	dest->size = 1;
	return;
}

void zClear(bignum * n)
{
	int i;
	for (i = 0; i <= NWORDS; i++)
		n->data[i] = 0;
	n->size = 1;
	return;
}

void zClearFull(bignum * n)
{
    int i;
    memset(n->data, 0, 2 * NWORDS * sizeof(base_t));
    n->size = 1;
    return;
}

bignum * zInit(void)
{
	int i;
	size_t sz = 2 * (NWORDS + 4);
	bignum *n;
	
	n = (bignum *)malloc(sizeof(bignum));

	n->data = (base_t *)xmalloc_align(sz * sizeof(base_t));
	for (i = 0; i < sz; i++)
	{
		n->data[i] = 0;
	}
	n->size = 1;

	return n;
}

void zFree(bignum *n)
{
	align_free(n->data);	
	free(n);
}

void zPrint(bignum *n)
{
    int j;
    for (j = MIN(n->size - 1, 2*NWORDS); j >= 0; j--)
        printf("%016lx", n->data[j]);
    return;
}

void zClamp(bignum * n)
{
	int j;
	int sn = abs(n->size);
    int sign = n->size < 0;

	for (j = sn - 1; j >= 0; j--)
	{
        if (n->data[j] == 0)
        {
            sn--;
        }
        else
            break;
	}

	n->size = (sn == 0 ? 1 : sn);	
	if (sign)
		n->size *= -1;

	return;
}

void zCopy(bignum * src, bignum * dest)
{
	//physically copy the digits of u into the digits of v
	int su = abs(src->size);
	int i;

	//memcpy(dest->data, src->data, su * sizeof(base_t));
	for (i = 0; i < su; i++)
	{
		dest->data[i] = src->data[i];
	}
	dest->size = src->size;
	return;
}

void zAdd(bignum * u, bignum * v, bignum * w)
{
	int i, su, sv;
	base_t *larger;
	base_t k;
	int n, m;

	if (u->size < 0)
	{
		if (v->size > 0)
		{
			//u is negative, v is not
			u->size *= -1;
			zSub(v, u, w);
			if (u != w)
				u->size *= -1;
			return;
		}
	}
	else if (v->size < 0)
	{
		//v is negative, u is not
		v->size *= -1;
		zSub(u, v, w);
		if (v != w)
			v->size *= -1;
		return;
	}

	su = abs(u->size);
	sv = abs(v->size);

	if (su >= sv)
	{
		larger = u->data;
		n = su;
		m = sv;
	}
	else
	{
		larger = v->data;
		n = sv;
		m = su;
	}

	k=0;
	for (i = 0; i < m; ++i)
		spAdd3(u->data[i], v->data[i], k, w->data + i, &k);

	for (     ; i < n; ++i)
		spAdd(larger[i], k, w->data + i, &k);

	w->size = n;
	if (k)
	{
		w->data[n] = k;
		w->size++;
	}

	// if one is negative then so is the other or we would be subtracting
	if (u->size < 0)
		w->size *= -1;

	return;
}

void zShortAdd(bignum * u, base_t v, bignum * w)
{
	int i, su;
	base_t k;

	if (u->size < 0)
	{
		//u is negative
		u->size *= -1;
		zShortSub(u, v, w);
		w->size *= -1;
		u->size *= -1;
		return;
	}

	su = abs(u->size);

	zCopy(u,w);

	//add
	spAdd(u->data[0], v, w->data, &k);

	//add the carry
	spAdd(u->data[1], k, w->data + 1, &k);

	if (k)
	{
		//only rarely will the carry propagate more than one place
		//special case this.
		for (i = 2; i < su; ++i)
			spAdd(u->data[i], k, w->data + i, &k);

		w->size = u->size;
		if (k)
		{
			w->data[u->size] = k;
			w->size++;
		}
	}

	return;
}

int zSub(bignum * u, bignum * v, bignum * w)
{
	base_t k = 0;
	int i, j, su, sv, sw, m, sign=0;
	base_t *bigger, *smaller;

	if (u->size < 0)
	{
		if (v->size > 0)
		{
			//u is negative, v is not, so really an addition
			u->size *= -1;
			zAdd(u, v, w);
			if (u != w)
				u->size *= -1;
			w->size *= -1;
            //printf("did an addition, result is neg\n");
			return 0;
		}
		else
		{
			//both are negative, so we really have -u + v or v - u
			v->size *= -1;
			u->size *= -1;
			zSub(v, u, w);
			if (v != w)
				v->size *= -1;
			if (u != w)
				u->size *= -1;
            //printf("both negative\n");
			return 0;
		}
	}
	else if (v->size < 0)
	{
		if (u->size > 0)
		{
			//v is negative, u is not, so really an addition
			v->size *= -1;
			zAdd(u, v, w);
			if (v != w)
				v->size *= -1;
            //printf("did an addition, result is pos\n");
			return 0;
		}
	}
	
	su = u->size;
	sv = v->size;

	if (su > sv) 
	{	
		bigger = u->data;
		smaller = v->data;
		sw = su;
		m = sv;
		goto beginsub;
	}
	if (su < sv)
	{
		bigger = v->data;
		smaller = u->data;
		sw = sv;
		m = su;
		sign=1;
		goto beginsub;
	}

	// same size
	m = su;
	sw = sv;
	for (i = su - 1; i >= 0; --i)
	{
		if (u->data[i] > v->data[i]) 
		{	
			bigger = u->data;
			smaller = v->data;
			goto beginsub;
		}
		if (u->data[i] < v->data[i])
		{	
			bigger = v->data;
			smaller = u->data;
			sign=1;
			goto beginsub;
		}
	}

	//equal if got to here
	w->size = 1;
	w->data[0] = 0;
	return 1;

beginsub:

	for (j = 0; j < m; ++j)
		spSub3(bigger[j], smaller[j], k, w->data + j, &k);


	//if there is a leftover word that is != 0, then subtract any
	//carry and simply copy any other leftover words

	//if there is a leftover word that is == 0, then subtract with
	//borrow for the rest of the leftover words.  this will happen rarely

	//leftover word?
	if (sw > m)
	{
		//not equal to zero?
		if (bigger[m] != 0)
		{
			//subtract any carry and copy the rest
			w->data[m] = bigger[m] - k;
			j = m + 1;
			m = sw;
			for(; j < m; j++)
				w->data[j] = bigger[j];
		}
		else
		{
			//equal to zero, need to subtract with borrow for the rest
			//of the leftover words.  
			j = m;
			m = sw;

			for ( ; j < m; ++j)
				spSub3(bigger[j], 0, k, w->data + j, &k);
		}
	}
	
	w->size = sw;
	zClamp(w);

	if (sign)
		w->size *= -1;

	if (w->size == 0)
		w->size = 1;

	return 0;	
}

void zShortSub(bignum * u, base_t v, bignum * w)
{
	// w = u - v
	// assume both are initially positive; result can be negative

	int i, su = abs(u->size);
	base_t k = 0;

	su = abs(u->size);
	w->size = su;

	if (u->size < 0)
	{
		//u is negative, really an addition
		u->size *= -1;
		zShortAdd(u,v,w);
		u->size *= -1;
		w->size *= -1;
		return;
	}

	zCopy(u,w);

	//subtract
	spSub3(u->data[0],v,0,w->data,&k);

	//subtract the borrow
	spSub3(u->data[1],k,0,w->data+1,&k);
	
	if (k)
	{
		//propagate the borrow
		for (i=2;i<su;++i)
			spSub3(u->data[i],0,k,w->data+i,&k);
	}

	//check if we lost the high digit
	if ((w->data[su - 1] == 0) && (su != 1))
		su--;
	w->size = su;

	//check for u < v
	if (k)
	{
		//then u < v, and result is negative
		w->data[0] = ~w->data[0];
		w->data[0]++;
		w->size *= -1;
	}

	return;
}

int zCompare(bignum * u, bignum * v)
{
	//return 1 if u > v, -1 if u < v, 0 if equal
	int i,j,su,sv;

	i = u->size < 0;
	j = v->size < 0;

	su = abs(u->size);
	sv = abs(v->size);
	
	if (i > j) 
	{
		//v pos, u neg
		//make sure both are not zero
		if ((u->data[0] == 0) && (su == 1) && (v->data[0] == 0) && (sv == 1))
			return 0;
		else
			return -1;	
	}
	if (j > i) 
	{
		//u pos, v neg
		//make sure both are not zero
		if ((u->data[0] == 0) && (su == 1) && (v->data[0] == 0) && (sv == 1))
			return 0;
		else
			return 1;	
	}	

	//check obvious
	if (j)
	{	//both are negative
		if (su > sv) return -1;
		if (su < sv) return 1;
	}
	else
	{	//both are positive
		if (su > sv) return 1;
		if (su < sv) return -1;
	}

	//if the numbers are both negative, then we'll need to switch the return value
	for (i = su - 1; i>=0; --i)
	{
		if (u->data[i] > v->data[i]) 
			return (1 - 2*j);
		if (u->data[i] < v->data[i])
			return (-1 + 2*j);
	}

	//equal if got to here
	return 0;
}

int zCompare1(bignum * u, base_t v)
{
    // return 1 if u > v, -1 if u < v, 0 if equal.
    // single digit v is assumed to be positive.
    if (u->size < 0)
    {
        return -1;
    }
    else if (u->size > 1)
    {
        return 1;
    }
    else if (u->data[0] > v)
    {
        return 1;
    }
    else if (u->data[0] < v)
    {
        return -1;
    }
    else
    {
        return 0;
    }
}

base_t zShortDiv(bignum * u, base_t v, bignum * q)
{
	// q = u/v
	// return the remainder

	int su = abs(u->size);
    int sign = u->size < 0 ? 1 : 0;
	int i;
	base_t rem = 0;

	q->size = su;

	i = su - 1;
	if (u->data[i] < v) 
	{
		rem = u->data[i];
		q->data[i--] = 0;
	}

	while (i >= 0)
	{
		base_t quot1;

#if DIGITBITS == 64
		__asm__ ("divq %4"
			: "=a"(quot1),"=d"(rem)
			: "1"(rem), "0"(u->data[i]), "r"(v) );
#else
		__asm__ ("divl %4"
			: "=a"(quot1),"=d"(rem)
			: "1"(rem), "0"(u->data[i]), "r"(v) );
#endif

		q->data[i] = quot1;
		i--;
	}

	//the quotient could be one limb smaller than the input
	if ((q->data[q->size - 1] == 0) && (q->size != 1))
		q->size--;

	if (sign)
		q->size *= -1;

	return rem;
}

void zDiv(bignum * u, bignum * v, bignum * q, bignum * r)
{
	/*
	q = u \ v
	r = u mod v
	u is overwritten

	schoolbook long division.  see knuth TAOCP, vol. 2
	*/

	base_t v1=0,v2=0,d=0,k,qhat,rhat,uj2,tt[2],pp[2];
	int i,j,m,su,sv;
	int s =0,cmp,sdd,sd;
	unsigned int shift;
	base_t bitmask;
	
	su = abs(u->size);
	sv = abs(v->size);
	m = su-sv;

	//v > u, so just set q = 0 and r = u
	if (su < sv)
	{	
		q->size = 1;
		zCopy(u,r);

		return;
	}

	if (sv == 1)
	{
		r->data[0] = zShortDiv(u, v->data[0], q);
		r->size = 1;
		s = (v->size < 0);
		if (s)
		{
			q->size *= -1;
			r->size *= -1;
		}
		return;
	}

	//u and v are the same length
	if (su == sv)
	{	
		cmp = zCompare(u,v);
		//v > u, as above
		if (cmp < 0)
		{	
			q->size = 1;
			zCopy(u,r);
			return;
		} 
		else if (cmp == 0)	//v == u, so set q = 1 and r = 0
		{	
			q->size = 1;
			q->data[0] = 1;
			r->size = 1;
			r->data[0] = 0;
			return;
		}
	}

	//normalize v by left shifting until the high bit of v is set (v1 >= floor(2^31))
	bitmask = HIBITMASK;
	for (shift = 0; shift < DIGITBITS; ++shift)
	{
		if (v->data[sv-1] & bitmask)
			break;
		bitmask >>= 1;
	}

	//normalize v by shifting left (x2) shift number of times
	//overflow should never occur to v during normalization
	zShiftLeft(v,v,shift);

	//left shift u the same amount - may get an overflow here
	zShiftLeft(u,u,shift);
	if (abs(u->size) == su)
	{	//no overflow - force extra digit
		if (u->size < 0)
			u->size--;
		else
			u->size++;
		u->data[su] = 0;
		su++;
	}
	else
		su++;

	//copy first two digits of v to local variables for quick access
	v1=v->data[sv-1];
	v2=v->data[sv-2];
	
	sdd=0;
	sd=0;
	//main loop
	for (j=0;j<=m;++j)
	{
		//calculate qhat
		tt[1] = u->data[su-j-1];		//first digit of normalized u
		tt[0] = u->data[su-j-2];		//second digit of normalized u
		if (tt[1] == v1)
			qhat = MAXDIGIT;
		else
			spDivide(&qhat, &rhat, tt, v1);

		//quick check if qhat is too big based on our initial guess involving 
		//the first two digits of u and v.
		uj2 = u->data[su-j-3];

		while (1)
		{
			spMultiply(qhat,v1,&pp[0],&pp[1]);
			shortSubtract(tt,pp,tt);
			if (tt[1]) break;
			tt[1] = tt[0]; tt[0] = uj2; 
			spMultiply(qhat,v2,&pp[0],&pp[1]);
			i = shortCompare(pp,tt);  //p = v2*qhat, t = (uj*b+uj1-qhat*v1)*b + uj2

			if (i == 1)
				qhat--;
			else
				break;
		}
	
		//keep track of the significant digits
		if (qhat > 0)
		{
			sdd = sdd + 1 + sd;
			sd = 0;
		}
		else if (sdd != 0)
			sd++;

		//multiply and subtract, in situ
		k=0;
		for (i=0;i<sv;++i)
		{
			s=(int)(su-j-sv+i-1);
			spMultiply(v->data[i],qhat,&pp[0],&pp[1]);
			spAdd(pp[0],k,&tt[0],&tt[1]);
			u->data[s] = u->data[s] - tt[0]; 
			//check if this result is negative, remember the borrow for the next digit
			if (u->data[s] > (u->data[s] + tt[0]))
				k = pp[1] + tt[1] + 1;
			else
				k = pp[1] + tt[1];
		}
		
		//if the final carry is bigger than the most significant digit of u, then qhat
		//was too big, i.e. qhat[v1v2...vn] > [u0u1u2...un]
		if (k > u->data[su-j-1])
		{
			//correct by decrementing qhat and adding back [v1v2...vn] to [u0u1...un]
			qhat--;
			//first subtract the final carry, yielding a negative number for [u0u1...un]
			u->data[su-j-1] -= k;
			//then add back v
			k=0;
			for (i=0;i<sv;i++)
				spAdd3(u->data[su-j-sv+i-1],v->data[i],k,&u->data[su-j-sv+i-1],&k);
			u->data[su-j-1] += k;
		}
		else //else qhat was ok, subtract the final carry
			u->data[su-j-1] -= k;

		//set digit of q
		q->data[m-j] = qhat;
	}
	q->size = sdd+sd;
	zCopy(u,r);

	for (s=r->size - 1; s>=0; --s)
	{
		if ((r->data[s] == 0) && (r->size > 0))
			r->size--;
		else
			break;
	}

	//unnormalize.
	zShiftRight(v,v,shift);
	zShiftRight(r,r,shift);

	s = (u->size < 0) ^ (v->size < 0);
	if (s)
	{
		q->size *= -1;
		r->size *= -1;
	}

	return;
}

int shortCompare(base_t p[2], base_t t[2])
{
	//utility function used in zDiv
	int i;

	for (i=1;i>=0;--i)
	{
		if (p[i] > t[i]) return 1;
		if (p[i] < t[i]) return -1;
	}
	return 0;
}

int shortSubtract(base_t u[2], base_t v[2], base_t w[2])
{
	//utility function used in zDiv
	base_t j=0;

	w[0] = u[0] - v[0];
	if (w[0] > (MAXDIGIT - v[0]))
	{
		j=1;
		w[0] = w[0] + MAXDIGIT + 1;
	}
	w[1] = u[1] - v[1] - j;
	
	return 1;
}

void zMult(bignum * u, bignum * v, bignum * w, bignum *tmp)
{
    //w = u*v
    base_t k = 0;
    int su, sv, i, j, signu, signv;
    base_t *wptr;
    int words = u->size;

    signu = u->size < 0;
    signv = v->size < 0;

    su = abs(u->size);
    sv = abs(v->size);

    //for each digit of u
    for (i = 0; i < su; ++i)
    {
        //take an inner product and add in-situ with the previous inner products
        k = 0;
        wptr = &tmp->data[i];
        for (j = 0; j < sv; ++j)
        {
            spMulAdd(u->data[i], v->data[j], wptr[j], k, &wptr[j], &k);
        }
        wptr[j] += k;
    }
    tmp->size = su + sv;

    zClamp(tmp);

    if (((u->size == 1) && (u->data[0] == 0)) || ((v->size == 1) && (v->data[0] == 0)))
    {
        w->size = 1;
        w->data[0] = 0;
    }
    else
    {
        zCopy(tmp, w);

        if (signu ^ signv)
            w->size *= -1;
    }

    return;
}

void zMul(bignum * u, bignum * v, bignum * w)
{
	//w = u*v
	base_t k = 0;
	int su, sv, i, j, signu, signv;
	base_t *wptr;
	int words = u->size;
	bignum *tmp;

	signu = u->size < 0;
	signv = v->size < 0;

	tmp = zInit();
	
	su = abs(u->size);
	sv = abs(v->size);

	//for each digit of u
	for (i = 0; i < su; ++i)
	{
		//take an inner product and add in-situ with the previous inner products
		k=0;
		wptr = &tmp->data[i];
		for (j = 0; j < sv; ++j)
		{
			spMulAdd(u->data[i], v->data[j], wptr[j], k, &wptr[j], &k);
		}
		wptr[j] += k;
	}
	tmp->size = su+sv;
	
	zClamp(tmp);

	if (((u->size == 1) && (u->data[0] == 0)) || ((v->size == 1) && (v->data[0] == 0)))
	{
		w->size = 1;
		w->data[0] = 0;
	}
	else
	{
		zCopy(tmp, w);

		if (signu ^ signv)
			w->size *= -1;
	}

	zFree(tmp);
	return;
}

void zModMul(bignum * u, bignum * v, bignum * n, bignum * w)
{
	bignum * t1, *t2;

	t1 = zInit();
	t2 = zInit();

	zMul(u,v,t1);
	zDiv(t1,n,t2,w);

	zFree(t1);
	zFree(t2);
	return;
}

void zModMuls(bignum * u, bignum * v, bignum * n, bignum * w, bignum *s1, bignum *s2)
{
    zMul(u, v, s1);
    zDiv(s1, n, s2, w);
    return;
}

void zModExp(bignum *d, bignum *b, bignum *e, bignum *m)
{
    // d = b^e mod m
    // all b and e vector elements can be different.
    // all m elements are the same.
    int i, word = 0, bit = 0;
    int j;

    bignum *s1, *s2, *bb, *t;

    s1 = zInit();
    s2 = zInit();
    bb = zInit();
    t = zInit();

    zCopy(b, bb);
    zSet1(d, 1);

    while (word < NWORDS)
    {
        if (e->data[word] & (1 << bit))
        {
            zModMuls(d, bb, m, d, s1, s2);
        }

        zModMuls(bb, bb, m, bb, s1, s2);

        bit++;
        if (bit == 32)
        {
            bit = 0;
            word++;
        }
    }

    zFree(s1);
    zFree(s2);
    zFree(bb);
    zFree(t);
    return;
}

void zShortMul(bignum * u, base_t v, bignum * w)
{
	//w = u * v
	//schoolbook multiplication, see knuth TAOCP, vol. 2
	base_t k=0;
	long i;
	long su;

	su = abs(u->size);

	//inner product
	for (i = 0; i < su; ++i)
		spMulAdd(u->data[i], v, 0, k, &w->data[i], &k);

	//if still have a carry, add a digit to w
	if (k)
	{
		w->data[su]=k;
		su++;
	}

	if (v == 0)
	{
		w->size = 1;
	}
	else
	{
		w->size = su;

		if (u->size < 0)
			w->size *= -1;
	}

	return;
}

void zSqr(bignum * x, bignum * w)
{
	//this routine is faster than the generic comba sqr on MSVC x86_32 builds.
	bignum *t;

	t = zInit();

	zCopy(x, t);
	zMul(x, x, w);

	zFree(t);

	return;
}

void zShiftLeft(bignum * a, bignum * b, int x)
{	
	/* Computes a = b << x */
	int i,wordshift;
	int y;
	int sb,j;
	base_t mask, carry, nextcarry;

	wordshift = x / DIGITBITS;
	x = x % DIGITBITS;

	//create a mask for the bits that will overflow each digit
	mask = HIBITMASK;
	for (i = 1; i < x; ++i)
		mask = (mask >> 1) | mask;

	if (x == 0) mask = 0x0;
	
	sb = abs(b->size);
	a->size = sb;

	//for each digit, remember the highest x bits using the mask, then shift.
	//the highest x bits becomes the lowest x bits for the next digit
	y = DIGITBITS - x;
	carry = 0;
	for (j = 0; j < sb; ++j)
	{
		nextcarry = (b->data[j] & mask) >> y;
		a->data[j] = (b->data[j] << x) | carry;
		carry = nextcarry;
	}
	
	if (carry)
	{
		a->data[sb] = carry;
		a->size++;
	}

	if (wordshift)
	{
		//now shift by any full words
		for (i=a->size - 1;i>=0;i--)
			a->data[i+wordshift] = a->data[i];
		//zero out the ones that were shifted
		for (i=wordshift-1;i>=0;i--)
			a->data[i] = 0;
		a->size += wordshift;
	}	

	if (b->size < 0)
		a->size *= -1;

	return;
}

void zShiftLeft_1(bignum * a, bignum * b)
{
    /* Computes a = b << 1 */
    int i;
    int y;
    int sb, j;
    base_t mask, carry, nextcarry;

    //create a mask for the bits that will overflow each digit
    mask = HIBITMASK;
    sb = abs(b->size);
    a->size = sb;

    //for each digit, remember the highest x bits using the mask, then shift.
    //the highest x bits becomes the lowest x bits for the next digit
    y = DIGITBITS - 1;
    carry = 0;
    for (j = 0; j < sb; ++j)
    {
        nextcarry = (b->data[j] & mask) >> y;
        a->data[j] = (b->data[j] << 1) | carry;
        carry = nextcarry;
    }

    if (carry)
    {
        a->data[sb] = carry;
        a->size++;
    }

    if (b->size < 0)
        a->size *= -1;

    return;
}

void zShiftRight(bignum * a, bignum * b, int x)
{	/* Computes a = b >> x */
	int i, y, sign, wordshift;
	int sb;
	base_t mask, carry, nextcarry;

	wordshift = x / DIGITBITS;
	x = x % DIGITBITS;

	//create a mask for the bits that will overflow each digit
	mask = 0x1;
	for (i = 1; i < x; ++i)
	{
		mask = (mask << 1) | mask;
	}
	if (x == 0) mask = 0x0;
	
    sign =( b->size < 0);
	sb = abs(b->size);
	a->size = sb;

	//for each digit, remember the lowest x bits using the mask, then shift.
	//the lowest x bits becomes the highest x bits for the next digit
	y = DIGITBITS - x;
	carry = 0;
	for (i = sb - 1; i >= 0; --i)
	{
		nextcarry = (b->data[i] & mask) << y;
		a->data[i] = b->data[i] >> x | carry;
		carry = nextcarry;
	}

	if ((a->data[sb-1] == 0) && (a->size > 1))
		a->size--;

	if (wordshift)
	{
		//now shift by any full words
		for (i=0;i<a->size - 1;i++)
			a->data[i] = a->data[i+wordshift];
		//zero out the ones that were shifted
		a->size -= wordshift;
	}	

	if (sign)
		a->size *= -1;

	return;
}

void zShiftRight_1(bignum * a, bignum * b)
{	/* Computes a = b >> x */
    int i, sign;
    int sb;
    base_t mask, carry, nextcarry;

    //create a mask for the bits that will overflow each digit
    mask = 0x1;

    sign = (b->size < 0);
    sb = abs(b->size);
    a->size = sb;

    //for each digit, remember the lowest x bits using the mask, then shift.
    //the lowest x bits becomes the highest x bits for the next digit
    carry = 0;
    for (i = sb - 1; i >= 0; --i)
    {
        nextcarry = (b->data[i] & mask) << 31;
        a->data[i] = (b->data[i] >> 1) | carry;
        carry = nextcarry;
    }

    if ((a->data[sb - 1] == 0) && (a->size > 1))
        a->size--;

    if (sign)
        a->size *= -1;

    return;
}

int zLEGCD(bignum *u, bignum *v, bignum *w)
{
	//use the Lehman-Euclid algorithm to calculate GCD(u,v) = w
	//Algorithm L in Knuth, 4.5.2 p. 329
	//assumes u,v nonnegative

	base_t aa,bb,cc,dd;
	int i,j,k,it;
	base_signed_t a,b,c,d,t;
	base_t up,vdp, q1, q2;
	base_t mask;
	bignum *y, *zz;		//t and w, in knuth
	bignum *x;			//tmp variable
	bignum *uu, *vv;	//so u and v don't get destroyed
	bignum *uh, *vh;
	base_t udp[2],vp[2];
    

#if DIGITBITS == 32
    mask = 0xff000000;
#else
	mask = 0xff00000000000000;
#endif

	i = zCompare1(u,0);
	j = zCompare1(v,0);

	if (i == 0) 
	{
		zCopy(v,w);
		return 1;
	}
	if (j == 0)
	{
		zCopy(u,w);
		return 1;
	}

	//temp variables should be twice as big as the input, to make room
	//for intermediate operations.  w should be as big as the biggest input.
	i = u->size;
	j = v->size;
	if (j > i)
		i = j;

	y = zInit();
	zz = zInit();
	x = zInit();
	uu = zInit();
	vv = zInit();
	uh = zInit();
	vh = zInit();

	//put bigger number in uu, other in vv.
	i = zCompare(u,v);
	if (i >= 0)
	{
		zCopy(u,uu);
		zCopy(v,vv);
	}
	else
	{
		zCopy(v,uu);
		zCopy(u,vv);
	}

	j=0;
	while (vv->size > 1)
	{
		//Step L1
		for (it=vv->size;it<uu->size;it++)
			vv->data[it]=0;
		vv->size = uu->size;
		//get the most significant 32 bits of u and v, such that uhat >= vhat
		uh->data[1] = uu->data[uu->size - 1];
		uh->data[0] = uu->data[uu->size - 2];
		vh->data[1] = vv->data[vv->size - 1];
		vh->data[0] = vv->data[vv->size - 2];
		uh->size = vh->size = 2;

		//rightshift until uhat is a single word
		//0xff000000 is magic
		if ((uh->data[1] & mask) > 0)
		{
			uh->data[0] = uh->data[1];
			vh->data[0] = vh->data[1];
		}
		else
		{
			i=0;
			aa=uh->data[1];
			while ((aa & MAXDIGIT) != 0)
			{
				aa >>= 1;
				i++;
			}
			zShiftRight(uh,uh,i);
			zShiftRight(vh,vh,i);
		}
		
		//make u',v',u'',v''
		up = uh->data[0];
		vdp = vh->data[0];
		if (up == MAXDIGIT)
		{
			udp[0] = 0;
			udp[1] = 1;
		}
		else
		{
			udp[0] = up+1;
			udp[1] = 0;
		}

		if (vdp == MAXDIGIT)
		{
			vp[0] = 0;
			vp[1] = 1;
		}
		else 
		{
			vp[0] = vdp+1;
			vp[1] = 0;
		}

		a=1; b=0; c=0; d=1; 
		
		k=0;
		while (1)
		{
			//Step L2:
			/*test quotient, protecting for overflow.  the conditions:
			0 <= uhat + a <= 2^32
			0 <= uhat + b < 2^32
			0 <= vhat + c < 2^32
			0 <= vhat + d <= 2^32
			will always hold.  hence only need to check for the case where
			uhat == MAX_DIGIT and a = 1 or vhat == MAX_DIGIT and d = 1
			*/

			//first check for /0
			if (((vp[0] == 0) && (vp[1] == 0)) || (vdp == 0))
				break;

			//u''/v''
			if (udp[1] == 1)
			{
				spDivide(&q2,(base_t *)&t,udp,vdp);
			}
			else
				q2 = udp[0]/vdp;

			//u'/v'
			if (vp[1] >= 1)
				q1=0;
			else
				q1 = up/vp[0];

			if (q1 != q2)
				break;

			//Step L3: Emulate Euclid
			t=a-q1*c;
			a=c;
			c=t;
			t=b-q1*d;
			b=d;
			d=t;
			t=up-q1*vp[0];
			up=vp[0]; 
			vp[0]=t;
			t=udp[0]-q2*vdp;
			udp[0]=vdp;
			vdp=t;
			k++;
			if (k>10000)
				goto free;
		}

		//Step L4: multiprecision step
		if (b==0)
		{
			for (i=vv->size-1;i>=0;i--)
			{
				if (vv->data[i] == 0)
					vv->size--;
				else
					break;
			}
			zDiv(uu,vv,y,x);	//y = u mod v
			zCopy(vv,uu);			//u = v
			zCopy(x,vv);			//v = y
		}
		else
		{
			//aa=abs(a);
			//bb=abs(b);
			//cc=abs(c);
			//dd=abs(d);
			if (a<0)
				aa = -a;
			else
				aa = a;
			if (b<0)
				bb = -b;
			else
				bb = b;
			if (c<0)
				cc = -c;
			else
				cc = c;
			if (d<0)
				dd = -d;
			else
				dd = d;
			zShortMul(uu,aa,y);	//y = A*u
			zShortMul(vv,bb,x);	//y = y + B*v
			if (a<0)
			{
				zSub(x,y,x);
				zCopy(x,y);
			}
			else if (b<0)
				zSub(y,x,y);
			else
				zAdd(y,x,y);				
			
			zShortMul(uu,cc,zz);	//z = c*u
			zShortMul(vv,dd,x);	//z = z + d*v
			if (c<0)
			{
				zSub(x,zz,x);
				zCopy(x,zz);
			}
			else // (d<0)
				zSub(zz,x,zz);

			zCopy(y,uu);				//u = y;
			zCopy(zz,vv);				//v = z;
		}
		j++;
		if (j>10000)
			goto free;
	}

	//here, the size of v is 1, so finish up with regular GCD
	zBinGCD(uu,vv,w);

free:
	zFree(y);
	zFree(zz);
	zFree(x);
	zFree(uu);
	zFree(vv);
	zFree(uh);
	zFree(vh);
	return 1;
}

int zBinGCD(bignum *u, bignum *v, bignum *w)
{
	//computes w = gcd(u,v) 
	//follows algorithm B.  p.321 Knuth Vol. 2

	bignum *uu, *vv, *t;
	long i=0,j;
	int k,sz;

	sz = abs(u->size);
	if (abs(v->size) > sz)
	{
		sz = abs(v->size);
	}

	i = zCompare1(u, 0);
	j = zCompare1(v, 0);
	if (i == 0) 
	{
		zCopy(v,w);
		return 1;
	}
	if (j == 0)
	{
		zCopy(u,w);
		return 1;
	}
		
	uu = zInit();
	vv = zInit();
	t = zInit();

	zCopy(u,uu);
	zCopy(v,vv);

	//find power of 2 such that u and v are not both even
	k = 0;
	while(((uu->data[0] & 0x1) == 0) && ((vv->data[0] & 0x1) == 0))
	{
		zShiftRight(uu,uu,1);
		zShiftRight(vv,vv,1);
		k++;
	}

	j=0;
	do
	{
		if ((uu->data[0] & 0x1) == 0)
			zShiftRight(uu,uu,1);
		else if ((vv->data[0] & 0x1) == 0)
			zShiftRight(vv,vv,1);
		else
		{
			zSub(uu,vv,t);
			zShiftRight(t,t,1);
			if (zCompare(uu,vv) < 0)
				zCopy(t,vv);
			else
				zCopy(t,uu);
		}
		++j;
		if (j>= 10000)
			break;
	} while (zCompare1(uu, 0) > 0);

	zClear(w);
	w->data[0] = 1;
	zShiftLeft(w,w,k);
	zMul(w,vv,uu);
	zCopy(uu,w);
	
	zFree(uu);
	zFree(vv);
	zFree(t);
	return j;
}

void xGCD(bignum *a, bignum *b, bignum *x, bignum *y, bignum *g)
{
	//compute the extended GCD of a, b, returning g = GCD(a,b) and x, y 
	//such that ax + by = GCD(a,b) if a,b are coprime
	bignum *t1, *t2, *t3, *u, *v, *r, *R, *q, *tmp;

//	int i;
	/*

	Step 1: 
	if a < b then 
	Set u=0, v=1, and r=b 
	Set U=1, V=0, and R=a 
	else 
	Set u=1, v=0, and r=a 
	Set U=0, V=1, and R=b 

	Step 2: 
	if R = 0 then return r (for the gcd) and no inverses exist. 
	if R = 1 then return R (for the gcd), V (for the inverse a(mod b)) and U (for the inverse of b(mod a)). 
	
	Step 3: 
	Calculate q = int(r/R) 
	Calculate t1 = u - U*q 
	Calculate t2 = v - V*q 
	Calculate t3 = r - R*q 
	set u=U, v=V, r=R 
	set U=t1, V=t2, R=t3 
	goto Step 2. 
	*/

	tmp = zInit();
	t1 = zInit();
	t2 = zInit();
	t3 = zInit();
	q = zInit();
	r = zInit();
	R = zInit();
	u = zInit();
	v = zInit();

	//need to check for temp allocation

	zClear(x);
	zClear(y);


	if (zCompare(a,b) < 0)
	{
		u->data[0]=0;
		v->data[0]=1;
		zCopy(b,r);
		x->data[0]=1;
		y->data[0]=0;
		zCopy(a,R);
	}
	else
	{
		u->data[0]=1;
		v->data[0]=0;
		zCopy(a,r);
		x->data[0]=0;
		y->data[0]=1;
		zCopy(b,R);
	}

	while (1)
	{
		if (zCompare1(R, 0) == 0)
		{
			zCopy(r,g);
            x->data[0] = 0;
            x->size = 1;
            y->data[0] = 0;
            y->size = 1;
			break;
		}

		if (zCompare1(R, 1) == 0)
		{
			zCopy(R,g);
			break;
		}

		zCopy(r,tmp);
		zDiv(tmp,R,q,t3);		//q = int(r/R), t3 = r % R

		zMul(q,x,tmp);			//t1 = u - U*q
		zSub(u,tmp,t1);

		zMul(q,y,tmp);			//t2 = v - V*q
		zSub(v,tmp,t2);

		zCopy(x,u);
		zCopy(y,v);
		zCopy(R,r);

		zCopy(t1,x);
		zCopy(t2,y);
		zCopy(t3,R);

		//printf("iteration %d: x = %s\n", i, z2decstr(x));
		//printf("iteration %d: y = %s\n", i, z2decstr(y));
		//printf("iteration %d: g = %s\n", i, z2decstr(g));
		//printf("iteration %d: r = %s\n", i, z2decstr(r));
		//printf("iteration %d: R = %s\n", i, z2decstr(R));
		//printf("iteration %d: q = %s\n", i, z2decstr(q));
		//printf("iteration %d: u = %s\n", i, z2decstr(u));
		//printf("iteration %d: v = %s\n", i, z2decstr(v));
		//i++;
	}

	if (x->size < 0)
	{
		x->size *= -1;
		zSub(b,x,x);
	}

	if (y->size < 0)
	{
		y->size *= -1;
		zSub(a,y,y);
	}

	zFree(tmp);
	zFree(t1);
	zFree(t2);
	zFree(t3);
	zFree(q);
	zFree(r);
	zFree(R);
	zFree(u);
	zFree(v);
	return;
}

void str2hexz(char in[], bignum * u)
{
	// convert a string to a bigint
	char *s2,*s;
	char **ptr = NULL;

    // assume input is base10, we convert 9 digits at a time (32 bit words)
	int i,j,su,base=10,step=9;
	bignum * t;

    // allocate space for a temporary bignum
	t = zInit();

	// work with a copy of in (because the first step in the conversion process
    // inserts null characters into the string...).  This could probably be changed.
	s = (char *)malloc(8192*sizeof(char));
	strcpy(s,in);

    // compute how many 9-digit decimal words we have in the string
	su = strlen(s)/step + (strlen(s)%step != 0);

	// read 9 characters of s at a time into a base-10 bignum, 'u'
	j=0;
	for (i=0;i<su-1;i++)
	{
		s2 = &s[strlen(s)] - step;
		t->data[j] = strtoul(s2,ptr,base);
		s2[0] = '\0';
		j++;
	}

	if (strlen(s) > 0)
	{
		s2 = s;
		ptr = &s2;
		t->data[j] = strtoul(s2,ptr,base);
	}
	t->size = j+1;

	// now convert the base-10 bignum to a binary bignum
	zDec2Hex(t,u);

    // clear the upper words, if any
    for (i = u->size; i < NWORDS; i++)
    {
        u->data[i] = 0;
    }

	free(s);
	zFree(t);
	return;
}

void zDec2Hex(bignum * u, bignum * v)
{
	// convert u[] in dec to v[] in hex by multiplying the ith digit by (1e9)*i
	// and adding to the previous digits

	bignum * a, *b, *vv;
	base_t d = MAX_DEC_WORD;
	int i, j;

	a = zInit();
	b = zInit();
	vv = zInit();
	zClear(v);

	a->data[0] = 1;
	for (i = 0; i < u->size; i++)
	{
		zShortMul(a, u->data[i], b);
		zAdd(vv, b, vv);
		zShortMul(a, d, a);
	}

	zClamp(vv);
	zCopy(vv, v);

	zFree(a);
	zFree(b);
	zFree(vv);

	return;
}

void zHex2Dec(bignum * u, bignum * v)
{
	//convert u[] in hex to v[] in decimal by repeatedly dividing
	//u by 1e9 = 0x3b9aca00
	//the remainder of the ith division is the ith decimal digit.
	//when the quotient = 0, stop

	bignum * a, *b;
	base_t d = MAX_DEC_WORD;
	base_t r = 0;
	int su = u->size;
	//because decimal takes more room than hex to store

	a = zInit();
	b = zInit();
	zClear(v);

	zCopy(u, a);
	v->size = 1;
	do
	{
		r = zShortDiv(a, d, b);
		v->data[v->size - 1] = r;
		v->size++;
		zCopy(b, a);
	} while (zCompare1(a, 0) != 0);
	v->size--;

	zFree(a);
	zFree(b);
	return;
}

char *z2decstr(bignum * n)
{
	//pass in a pointer to a string.  if necessary, this routine will 
	//reallocate space for the string to accomodate its size.  If this happens
	//the pointer to the string's (likely) new location is automatically
	//updated and returned.
	bignum * a;
	int i,sza,sign = 0;
	char *tmp, *s;
	int nchars, j;

	a = zInit();

	s = (char *)malloc(8192 * sizeof(char));

	strcpy(s,"");
    nchars = 1;
    if (n->size < 0)
    {
        sign = 1;
        n->size *= -1;
        sprintf(s, "-");
        nchars++;
    }

	zHex2Dec(n, a);
	sza = abs(a->size);
	
	tmp = (char *)malloc((DEC_DIGIT_PER_WORD + 10) * sizeof(char));

	//print first word
#if DIGITBITS == 64
	sprintf(s,"%s%lu", s, a->data[sza - 1]);
#else
    sprintf(s, "%s%u", s, a->data[sza - 1]);
#endif
	nchars += ndigits_1(a->data[sza-1]) - 1;

	//print the rest
	for (i=sza - 2; i >= 0; i--)
	{
#if DIGITBITS == 64
		sprintf(tmp,"%019lu",a->data[i]);
#else
        sprintf(tmp, "%09u", a->data[i]);
#endif
		memcpy(s + nchars, tmp, DEC_DIGIT_PER_WORD * sizeof(char));
		nchars += DEC_DIGIT_PER_WORD;
	}
	s[nchars] = '\0';

	free(tmp);
	zFree(a);

    if (sign)
    {
        n->size *= -1;
    }

	return s;
}

void spMulAdd(base_t u, base_t v, base_t w, base_t t, base_t *lower, base_t *carry)
{
	base_t k,p;
	spMultiply(u,v,&p,carry);
	spAdd3(p,w,t,lower,&k);
	*carry += k;
	return;
}

void spMulMod(base_t u, base_t v, base_t m, base_t *w)
{
	base_t p[2];
	base_t q;

	spMultiply(u,v,&p[0],&p[1]);
	spDivide(&q,w,p,m);

	return;
}

