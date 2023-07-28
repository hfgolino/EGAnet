/*  Written in 2019 by David Blackman and Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

#include <time.h>
#include <stdint.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include "nanotime.h"

/* This is xoshiro128++ 1.0, one of our 32-bit all-purpose, rock-solid
   generators. It has excellent speed, a state size (128 bits) that is
   large enough for mild parallelism, and it passes all tests we are aware
   of.

   For generating just single-precision (i.e., 32-bit) floating-point
   numbers, xoshiro128+ is even faster.

   The state must be seeded so that it is not everywhere zero. */


static inline uint32_t rotl(const uint32_t x, int k) {
	return (x << k) | (x >> (32 - k));
}


static uint32_t s[4];

uint32_t next(void) {
	const uint32_t result = rotl(s[0] + s[3], 7) + s[0];

	const uint32_t t = s[1] << 9;

	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];

	s[2] ^= t;

	s[3] = rotl(s[3], 11);

	return result;
}


/* This is the jump function for the generator. It is equivalent
   to 2^64 calls to next(); it can be used to generate 2^64
   non-overlapping subsequences for parallel computations. */

void jump(void) {
	static const uint32_t JUMP[] = { 0x8764000b, 0xf542d2d3, 0x6fa035c3, 0x77f2db5b };

	uint32_t s0 = 0;
	uint32_t s1 = 0;
	uint32_t s2 = 0;
	uint32_t s3 = 0;
	for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
		for(int b = 0; b < 32; b++) {
			if (JUMP[i] & UINT32_C(1) << b) {
				s0 ^= s[0];
				s1 ^= s[1];
				s2 ^= s[2];
				s3 ^= s[3];
			}
			next();
		}

	s[0] = s0;
	s[1] = s1;
	s[2] = s2;
	s[3] = s3;
}


/* This is the long-jump function for the generator. It is equivalent to
   2^96 calls to next(); it can be used to generate 2^32 starting points,
   from each of which jump() will generate 2^32 non-overlapping
   subsequences for parallel distributed computations. */

void long_jump(void) {
	static const uint32_t LONG_JUMP[] = { 0xb523952e, 0x0b6f099f, 0xccf5a0ef, 0x1c580662 };

	uint32_t s0 = 0;
	uint32_t s1 = 0;
	uint32_t s2 = 0;
	uint32_t s3 = 0;
	for(int i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
		for(int b = 0; b < 32; b++) {
			if (LONG_JUMP[i] & UINT32_C(1) << b) {
				s0 ^= s[0];
				s1 ^= s[1];
				s2 ^= s[2];
				s3 ^= s[3];
			}
			next();
		}

	s[0] = s0;
	s[1] = s1;
	s[2] = s2;
	s[3] = s3;
}

// Function to get 4 random seeds
uint32_t splitmix64(uint32_t *x) {
    uint32_t z = (*x += 0x9e3779b97f4a7c15);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    return z ^ (z >> 31);
}

// Function to set single seed to get the 4 random seeds
void seed_xoshiro256(uint32_t seed) {
    for(int i = 0; i < 4; i++) {
        seed = splitmix64(&seed);
        s[i] = seed;
    }
}

// Function to generate random uniform data between 0 and 1
float xoshiro_uniform(void) {
  return (float) next() / ((float) UINT32_MAX + 1);
}


// Function get "random" seeds
SEXP r_xoshiro_seeds(SEXP n, SEXP r_seed) {

    // Initialize R values
    int n_values = INTEGER(n)[0];
    uint32_t seed_value = (uint32_t) REAL(r_seed)[0];

    // For random seed, use zero
    if(seed_value == 0) { // Use clocktime in nanoseconds
      seed_value = (uint32_t) get_time_ns();
    }

    // Seed the random number generator
    seed_xoshiro256(seed_value);

    // Create R vector
    SEXP r_output = PROTECT(allocVector(REALSXP, n_values));

    // Generate a random number and store it in the array
    for(int i = 0; i < n_values; i++) {
        REAL(r_output)[i] = (double) next();
    }

    // Release protected SEXP objects
    UNPROTECT(1);

    // Return the result
    return r_output;

}

// Function to shuffle indices without replacement (Fisher-Yates or Knuth shuffle)
SEXP r_xoshiro_shuffle(SEXP r_vector, SEXP r_seed) {

    // Get length of R vector
    int vector_length = length(r_vector);

    // Initialize R values
    uint32_t seed_value = (uint32_t) REAL(r_seed)[0];

    // For random seed, use zero
    if(seed_value == 0) { // Use clocktime in nanoseconds
      seed_value = (uint32_t) get_time_ns();
    }

    // Seed the random number generator
    seed_xoshiro256(seed_value);

    // Protect the input SEXP
    PROTECT(r_vector);

    // Shuffle the array using the Fisher-Yates algorithm
    for (int i = vector_length - 1; i > 0; i--) {
        int j = next() % (i + 1); // generates random index between 0 and i
        int tmp = INTEGER(r_vector)[j];
        INTEGER(r_vector)[j] = INTEGER(r_vector)[i];
        INTEGER(r_vector)[i] = tmp;
    }

    // Unprotect the SEXP
    UNPROTECT(1);

    // Return the result
    return r_vector;

}

// Function to shuffle indices with replacement
SEXP r_xoshiro_shuffle_replace(SEXP r_vector, SEXP r_seed) {

    // Get length of R vector
    uint32_t vector_length = (uint32_t) length(r_vector);

    // Initialize R values
    uint32_t seed_value = (uint32_t) REAL(r_seed)[0];

    // For random seed, use zero
    if(seed_value == 0) { // Use clocktime in nanoseconds
      seed_value = (uint32_t) get_time_ns();
    }

    // Seed the random number generator
    seed_xoshiro256(seed_value);

    // Create R vector
    SEXP r_output = PROTECT(allocVector(REALSXP, vector_length));

    // Shuffle
    for(int i = 0; i < vector_length; i++) {
        REAL(r_output)[i] = (double) (next() % (vector_length)) + 1;
    }

    // Release protected SEXP objects
    UNPROTECT(1);

    // Return the result
    return r_output;

}
