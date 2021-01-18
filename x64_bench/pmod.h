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

#ifndef _PMOD_H
#define _PMOD_H

// Modular exponentiation relies on a big-integer math library and libraries
// that perform modular arithmetic.  We define routines that use a
// homegrown bigint library.
#include "bigarith.h"
#include "monty_arith.h"

typedef struct
{
    bignum **libpmod_gwin;
} pmod_t;

#define MAX_WINSIZE 8

int get_winsize(void);
int get_bitwin(bignum *b, int bitloc, int winsize, int winmask);
void lr_powm(pmod_t *pmod_state, monty *mdata, bignum *c, bignum *a, bignum *b, bignum *n, bignum *s);
void lrwin_powm(pmod_t *pmod_state, monty *mdata, bignum *c, bignum *a, bignum *b, bignum *n, bignum *s);
void lroddwin_powm(pmod_t *pmod_state, monty *mdata, bignum *c, bignum *a, bignum *b, bignum *n, bignum *s);

void pmodlib_init(pmod_t *pmod_state);
void pmodlib_free(pmod_t *pmod_state);



#endif