# 
# Copyright (c) 2014, Ben Buhrow
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer. 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# The views and conclusions contained in the software and documentation are those
# of the authors and should not be interpreted as representing official policies, 
# either expressed or implied, of the FreeBSD Project.
# 
# 
# Copyright (c) 2018 by The Mayo Clinic, though its Special Purpose
#  Processor Development Group (SPPDG). All Rights Reserved Worldwide.
#  Licensed under the Apache License, Version 2.0 (the "License"); you may
#  not use this file except in compliance with the License. You may obtain
#  a copy of the License at http://www.apache.org/licenses/LICENSE-2.0.
# Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied,
#            including conditions of title, non-infringement, merchantability,
#           or fitness for a particular purpose
#  See the License for the specific language governing permissions and
#  limitations under the License.
# This file is a snapshot of a work in progress, originated by Mayo
#  Clinic SPPDG.




#--------------------------- flags -------------------------
CC = icc
#CFLAGS = -g -march=core2 -mtune=core2
#CFLAGS = -static
#CFLAGS = -S -fsource-asm
WARN_FLAGS = -Wall #-W -Wconversion
OPT_FLAGS = -O2
INC = -I. 
LIBS =
BINNAME = avx512_modexp
CFLAGS += -I../gmp_install/gmp-6.2.0/include/ 
CFLAGS += -L../gmp_install/gmp-6.2.0/lib/ 
CFLAGS += -g -gdwarf-4

#--------------------------- make options -------------------------


ifeq ($(COMPILER),mingw)
# NOTE: Using -fcall-used instead of -ffixed is much better and still works.
# -fcall-used simply prevents the named registers from being saved/restored while
# -ffixed prevents them from being used at all.  The code benefits a lot from being
# able to use all 32 zmm registers.
	CC = gcc
    BINNAME = avx512_modexp_mingw
	CFLAGS += -fopenmp
    CFLAGS += -fcall-used-xmm16 -fcall-used-xmm17 -fcall-used-xmm18 -fcall-used-xmm19
    CFLAGS += -fcall-used-xmm20 -fcall-used-xmm21 -fcall-used-xmm22 -fcall-used-xmm23
    CFLAGS += -fcall-used-xmm24 -fcall-used-xmm25 -fcall-used-xmm26 -fcall-used-xmm27
    CFLAGS += -fcall-used-xmm28 -fcall-used-xmm29 -fcall-used-xmm30 -fcall-used-xmm31
else ifeq ($(COMPILER),gcc730)
    CC = gcc-7.3.0
	CFLAGS += -fopenmp 
else ifeq ($(COMPILER),gcc11)
    CC = gcc-11.1.0
	CFLAGS += -fopenmp 
else
	CFLAGS += -qopenmp
endif

ifdef MAXBITS
	CFLAGS += -DMAXBITS=$(MAXBITS)
endif

ifdef BASE52
	CFLAGS += -DBASE52
endif

ifeq ($(KNL),1)
    ifeq ($(COMPILER),gcc)
        CFLAGS += -march=knl -DTARGET_KNL
    else
        CFLAGS += -xMIC-AVX512 -DTARGET_KNL
    endif
    OBJ_EXT = .o
    BINNAME := ${BINNAME:%=%_knl}
else
    OBJ_EXT = .o
    ifeq ($(SKYLAKEX),1)
        CFLAGS += -DSKYLAKEX
        ifeq ($(COMPILER),icc)
            CFLAGS += -march=skylake-avx512  -DTARGET_KNL
        else
            CFLAGS += -march=skylake-avx512  -DTARGET_KNL
        endif
    else
        OPT_FLAGS += -mavx
    endif
endif


ifeq ($(CC),icc)
	CFLAGS += -qmkl
endif


ifeq ($(PROFILE),1)
	CFLAGS += -pg
	BINNAME := ${BINNAME:%=%_prof}
endif


CFLAGS += -g $(OPT_FLAGS) $(WARN_FLAGS) $(INC)

ifeq ($(STATIC),1)
	CFLAGS += -static-intel -static
	LIBS += -L/usr/lib/x86_64-redhat-linux6E/lib64/ -lm
else
	LIBS += -lm -lgmp
endif

	
#--------------------------- file lists -------------------------
SRCS = \
	common.c \
	vecarith52.c \
	vecarith.c \
	main.c 

OBJS = $(SRCS:.c=$(OBJ_EXT))



#---------------------------Header file lists -------------------------
HEAD = \
	vecarith.h

#---------------------------Make Targets -------------------------

all: $(OBJS)
	rm -f libvecarith.a
	ar r libvecarith.a $(OBJS)
	ranlib libvecarith.a
	$(CC) $(CFLAGS) $(OBJS) -o $(BINNAME) libvecarith.a $(LIBS)


clean:
	rm -f $(OBJS)
	
#---------------------------Build Rules -------------------------

	
%$(OBJ_EXT): %.c $(HEAD)
	$(CC) $(CFLAGS) -c -o $@ $<

