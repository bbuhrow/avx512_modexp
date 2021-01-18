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

#include "util.h"


// ============================================================================
// precision time
// ============================================================================
#ifdef _MSC_VER

	/* Core aware timing on Windows, courtesy of Brian Gladman */

	#if defined( _WIN64 )

		#define current_processor_number GetCurrentProcessorNumber

	#else

		unsigned long current_processor_number(void)
		{
			__asm
			{
				mov     eax,1
				cpuid
				shr     ebx,24
				mov     eax, ebx
			}
		}

	#endif

	int lock_thread_to_core(void)
	{   DWORD_PTR afp, afs;

		if(GetProcessAffinityMask(GetCurrentProcess(), &afp, &afs))
		{
			afp &= (DWORD_PTR)(1 << current_processor_number());
			if(SetThreadAffinityMask(GetCurrentThread(), afp))
				return EXIT_SUCCESS;
		}
		return EXIT_FAILURE;
	}

	int unlock_thread_from_core(void)
	{   DWORD_PTR afp, afs;

        if(GetProcessAffinityMask(GetCurrentProcess(), &afp, &afs))
		{
			if(SetThreadAffinityMask(GetCurrentThread(), afp))
				return EXIT_SUCCESS;
		}
		return EXIT_FAILURE;
	}

    double cycles_per_second = 0.0;
    double ticks_per_second = 0.0;
    double cycles_per_tick = 0.0;

    uint64_t measure_processor_speed(int millisec)
    {   unsigned long long cycles;

        lock_thread_to_core();
        cycles = __rdtsc();
        Sleep(millisec);
        cycles = __rdtsc() - cycles;
        unlock_thread_from_core();
        cycles_per_second = 10.0 * (double)cycles;

        if(ticks_per_second == 0.0)
        {   LARGE_INTEGER ll;
            QueryPerformanceFrequency(&ll);
            ticks_per_second = (double)ll.QuadPart;
            cycles_per_tick = cycles_per_second / ticks_per_second;
        }
        return cycles;
    }

    double get_tsc_time(void)
    {
        if(cycles_per_second == 0.0)
            measure_processor_speed(100); 
        return __rdtsc() / cycles_per_second;
    }

    double get_pfc_time(void)
    {   LARGE_INTEGER ll;

        if(ticks_per_second == 0.0)
            measure_processor_speed(100); 
        QueryPerformanceCounter(&ll);
        return ll.QuadPart / ticks_per_second;
    }

#else

	double cycles_per_second = 0.0;

	uint64_t measure_processor_speed(int millisec)
	{   
		uint64_t cycles;
		struct timeval start, stop;
		double t_time;

		gettimeofday(&start,NULL);

		cycles = read_clock(); 
		do
		{
			gettimeofday (&stop, NULL);
            t_time = my_difftime (&start, &stop);
		}
		while (t_time*1000 < millisec);
		cycles = read_clock() - cycles;

		return cycles;                  /* return cycles per second  */
	}

#endif


uint64_t read_clock(void) 
{
#if defined(__GNUC__) && (defined(__i386__) || defined(GCC_ASM64X) )
	uint32_t lo, hi;
	asm("rdtsc":"=d"(hi),"=a"(lo));
	return (uint64)hi << 32 | lo;

#elif defined(_MSC_VER)
    LARGE_INTEGER ll;
    QueryPerformanceCounter(&ll);
	return (uint64_t)(ll.QuadPart * cycles_per_tick);
#else
	struct timeval thistime;   
	gettimeofday(&thistime, NULL);
	return (uint64_t)(cycles_per_second * 
        (thistime.tv_sec + thistime.tv_usec / 1000000.0));
#endif
}

#ifdef _MSC_VER
int gettimeofday(struct timeval *tv, struct timezone *tz)
{
	FILETIME ft;
	unsigned __int64 tmpres = 0;
	static int tzflag;
	 
	if (NULL != tv)
	{
	GetSystemTimeAsFileTime(&ft);
	 
	tmpres |= ft.dwHighDateTime;
	tmpres <<= 32;
	tmpres |= ft.dwLowDateTime;
	 
	/*converting file time to unix epoch*/
	tmpres /= 10;  /*convert into microseconds*/
	tmpres -= DELTA_EPOCH_IN_MICROSECS; 
	tv->tv_sec = (long)(tmpres / 1000000UL);
	tv->tv_usec = (long)(tmpres % 1000000UL);
	}
	 
	if (NULL != tz)
	{
	if (!tzflag)
	{
		_tzset();
		tzflag++;
	}
	tz->tz_minuteswest = _timezone / 60;
	tz->tz_dsttime = _daylight;
	}
	 
	return 0;
}
#endif

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


// ============================================================================
// randomness
// ============================================================================
void get_random_seeds(rand_t *r) {

	uint32_t tmp_seed1, tmp_seed2;

#ifndef WIN32

	FILE *rand_device = fopen("/dev/urandom", "r");

	if (rand_device != NULL) {
		fread(&tmp_seed1, sizeof(uint32_t), (size_t)1, rand_device);
		fread(&tmp_seed2, sizeof(uint32_t), (size_t)1, rand_device);
		fclose(rand_device);
	}
	else

#endif
	{
		/* <Shrug> For everyone else, sample the current time,
		   the high-res timer (hopefully not correlated to the
		   current time), and the process ID. Multithreaded
		   applications should fold in the thread ID too */

		uint64_t high_res_time = read_clock();
		tmp_seed1 = ((uint32_t)(high_res_time >> 32) ^
			     (uint32_t)time(NULL)) * 
			    (uint32_t)getpid();
		tmp_seed2 = (uint32_t)high_res_time;
	}

	/* The final seeds are the result of a multiplicative
	   hash of the initial seeds */

	r->low = tmp_seed1 * ((uint32_t)40499 * 65543);
	r->hi = tmp_seed2 * ((uint32_t)40499 * 65543);
}

// Knuth's 64 bit MMIX LCG, using a global 64 bit state variable.
uint32_t spRand(uint64_t *state, uint32_t lower, uint32_t upper)
{
	// advance the state of the LCG and return the appropriate result
    *state = 6364136223846793005ULL * (*state) + 1442695040888963407ULL;
	return lower + (uint32_t)(
        (double)(upper - lower) * (double)(*state >> 32) * INV_2_POW_32);
}

uint64_t spRand64_range(uint64_t *state, uint64_t lower, uint64_t upper)
{
    // advance the state of the LCG and return the appropriate result
    *state = 6364136223846793005ULL * (*state) + 1442695040888963407ULL;
    return lower + (uint64_t)(
        (double)(upper - lower) *  ((double)*state * INV_2_POW_64));
}

uint64_t spRand64(uint64_t *state)
{
    // advance the state of the LCG and return the appropriate result.
    // assume lower = 0 and upper = maxint
    *state = 6364136223846793005ULL * (*state) + 1442695040888963407ULL;
    return *state;
}

// ============================================================================
// hashing
// ============================================================================

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

