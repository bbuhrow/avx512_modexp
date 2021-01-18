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


#--------------------------- EXAMPLE BUILDS -------------------------
# make MAXBITS=256 
# make MAXBITS=512 



#--------------------------- flags -------------------------
CC = gcc
WARN_FLAGS = -Wall #-W -Wconversion
OPT_FLAGS = -O3
INC = -I. -I../gmp-6.1.2/include/
LIBS = -L../gmp-6.1.2/lib/
BINNAME = pmod_bench

#--------------------------- make options -------------------------


ifeq ($(COMPILER),gcc73)
	CC = gcc-7.3.0
endif

ifdef MAXBITS
	CFLAGS += -DMAXBITS=$(MAXBITS)
endif

ifdef VERBOSE
  CFLAGS += -DVERBOSE=$(VERBOSE)
endif

CFLAGS += -g $(OPT_FLAGS) $(WARN_FLAGS) $(INC)
LIBS += -lm -lgmp
	
#--------------------------- file lists -------------------------
SRCS = \
	bigarith.c \
	x64_arith.c \
	monty_arith.c \
	pmod.c \
	main.c \
	util.c

OBJS = $(SRCS:.c=.o)

#---------------------------Header file lists -------------------------
HEAD = \
	monty_arith.h \
	bigarith.h \
	x64_arith.h \
	pmod.h \
	util.h

#---------------------------Make Targets -------------------------

all: $(OBJS)
	rm -f libpmod.a
	ar r libpmod.a $(OBJS)
	ranlib libpmod.a
	$(CC) $(CFLAGS) $(OBJS) -o $(BINNAME) libpmod.a $(LIBS)


clean:
	rm -f $(OBJS)
	
#---------------------------Build Rules -------------------------

	
%$(OBJ_EXT): %.c $(HEAD)
	$(CC) $(CFLAGS) -c -o $@ $<

