/*  Written in 2019 by David Blackman and Sebastiano Vigna (vigna@acm.org)

 To the extent possible under law, the author has dedicated all copyright
 and related and neighboring rights to this software to the public domain
 worldwide. This software is distributed without any warranty.

 See <http://creativecommons.org/publicdomain/zero/1.0/>. */

#include <stdint.h>
#include <R.h>
#include <Rinternals.h>
#include "nanotime.h"

/* This is xoshiro256++ 1.0, one of our all-purpose, rock-solid generators.
 It has excellent (sub-ns) speed, a state (256 bits) that is large
 enough for any parallel application, and it passes all tests we are
 aware of.

 For generating just floating-point numbers, xoshiro256+ is even faster.

 The state must be seeded so that it is not everywhere zero. If you have
 a 64-bit seed, we suggest to seed a splitmix64 generator and use its
 output to fill s. */

static inline uint64_t rotl(const uint64_t x, int k) {
  return (x << k) | (x >> (64 - k));
}

static uint64_t s[4];

uint64_t next(void) {
  const uint64_t result = rotl(s[0] + s[3], 23) + s[0];

  const uint64_t t = s[1] << 17;

  s[2] ^= s[0];
  s[3] ^= s[1];
  s[1] ^= s[2];
  s[0] ^= s[3];

  s[2] ^= t;

  s[3] = rotl(s[3], 45);

  return result;
}

/* `jump` and `long_jump` have been removed from the original
    source file because they are not used in {EGAnet} */

// Get initial states
uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9e3779b97f4a7c15);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    return z ^ (z >> 31);
}

// Function to set single seed to get the 4 random seeds
void seed_xoshiro256(uint64_t seed) {
    for(int i = 0; i < 4; i++) {
        seed = splitmix64(&seed);
        s[i] = seed;
    }
}

// Function to generate random uniform data between 0 and 1
double xoshiro_uniform(void) {
  return (double) (next() / ((double) UINT64_MAX + 1));
}

// Function to uniform values into R
SEXP r_xoshiro_uniform(SEXP n, SEXP r_seed) {

  // Initialize R values
  int n_values = INTEGER(n)[0];
  uint64_t seed_value = (uint64_t) REAL(r_seed)[0];

  // For random seed, use zero
  if(seed_value == 0) { // Use clocktime in nanoseconds
    seed_value = get_time_ns();
  }

  // Seed the random number generator
  seed_xoshiro256(seed_value);

  // Create R vector
  SEXP r_output = PROTECT(allocVector(REALSXP, n_values));

  // Generate a random number and store it in the array
  for(int i = 0; i < n_values; i++) {
    REAL(r_output)[i] = xoshiro_uniform();
  }

  // Release protected SEXP objects
  UNPROTECT(1);

  // Return the result
  return r_output;

}

// Function to pre-generate seeds
SEXP r_xoshiro_seeds(SEXP n, SEXP r_seed) {

    // Initialize R values
    int n_values = INTEGER(n)[0];
    uint64_t seed_value = (uint64_t) REAL(r_seed)[0];

    // For random seed, use zero
    if(seed_value == 0) { // Use clocktime in nanoseconds
      seed_value = get_time_ns();
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
    uint64_t seed_value = (uint64_t) REAL(r_seed)[0];

    // For random seed, use zero
    if(seed_value == 0) { // Use clocktime in nanoseconds
      seed_value = get_time_ns();
    }

    // Seed the random number generator
    seed_xoshiro256(seed_value);

    // Protect the input SEXP
    PROTECT(r_vector);

    // Shuffle the array using the Fisher-Yates (or Knuth shuffle) algorithm
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
    int vector_length = length(r_vector);

    // Initialize R values
    uint64_t seed_value = (uint64_t) REAL(r_seed)[0];

    // For random seed, use zero
    if(seed_value == 0) { // Use clocktime in nanoseconds
      seed_value = get_time_ns();
    }

    // Seed the random number generator
    seed_xoshiro256(seed_value);

    // Create R vector
    SEXP r_output = PROTECT(allocVector(REALSXP, vector_length));

    // Shuffle
    for(int i = 0; i < vector_length; i++) {
        REAL(r_output)[i] = (double) (next() % vector_length) + 1;
    }

    // Release protected SEXP objects
    UNPROTECT(1);

    // Return the result
    return r_output;

}
