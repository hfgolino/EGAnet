/*  Written in 2019 by David Blackman and Sebastiano Vigna (vigna@acm.org)

 To the extent possible under law, the author has dedicated all copyright
 and related and neighboring rights to this software to the public domain
 worldwide. This software is distributed without any warranty.

 See <http://creativecommons.org/publicdomain/zero/1.0/>. */

#include <stdint.h>
#include <R.h>
#include <Rinternals.h>
#include "nanotime.h"
#include "xoshiro.h"

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

uint64_t next(xoshiro256_state* state) {
    const uint64_t result = rotl(state->s[0] + state->s[3], 23) + state->s[0];

    const uint64_t t = state->s[1] << 17;

    state->s[2] ^= state->s[0];
    state->s[3] ^= state->s[1];
    state->s[1] ^= state->s[2];
    state->s[0] ^= state->s[3];

    state->s[2] ^= t;

    state->s[3] = rotl(state->s[3], 45);

    return result;
}

/*

 `jump` and `long_jump` have been removed from the original
  source file because they are not used in {EGAnet}

 `jump` ensures non-overlapping subsequences within xoshiro256++;
 however, {EGAnet} does not take advantage of it

 instead, seeds are pre-generated and allow splitmix64 to
 do some of the work

 according to Blackman and Vigna (2019), splitmix64 passes BigCrush,
 meaning there is sufficient randomness when setting xoshiro256++'s state

 further, in a footnote (p.11), they state that "with 256 bits of state,
 2^64 sequences of length 2^64 starting at 2^64 random points in the
 state space have an overlap probability less than 2^-64"

 leaving the sequence length aside, there is never a case that people
 should ever get close to 2^64 random starting points (seeds) in {EGAnet}
 (e.g., 500 is the default for `bootEGA`) -- even with Monte Carlo
 simulation work most values wouldn't cross 2^20 or 1 million
 random starting points

 based on Vigna (2020), the upper bound on the overlap probability
 is formally defined as n^2 * L / P where n = processors (or starting points),
 L is sequence length, and P is the period of the PRNG

 using a sequence length of 2^64 (far beyond anything used in {EGAnet}
 and Monte Carlo simulations in quantitative psychology) and xoshiro256++'s
 period of 2^256 - 1, then some upper bound probabilities of n
 random starting points are defined below:

 500 (default `bootEGA`) = 3.98273e-53

 1e06 (large simulation) = 1.593092e-46

 1e12 (one trillion seeds) = 1.593092e-34

 taken together, `jump` ensures non-overlapping subsequences, but
 the approach applied here is common and the state space is large
 enough that although there is a non-zero chance of subsequence
 overlap, the result is extremely unlikely < 1.593092e-34

*/

// Get initial states
uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9e3779b97f4a7c15);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    return z ^ (z >> 31);
}

// Function to set single seed to get the 4 random seeds
void seed_xoshiro256(xoshiro256_state* state, uint64_t seed) {
    
    // Loop unrolled
    seed = splitmix64(&seed); state->s[0] = seed;
    seed = splitmix64(&seed); state->s[1] = seed;
    seed = splitmix64(&seed); state->s[2] = seed;
    seed = splitmix64(&seed); state->s[3] = seed;
    
}

// Stored so it doesn't need to be repeatedly called
static const double UINT64_MAX_PLUS_ONE = (double) UINT64_MAX + 1;

// Function to generate random uniform data between 0 and 1
double xoshiro_uniform(xoshiro256_state* state) {
  return (double) (next(state) / UINT64_MAX_PLUS_ONE);
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
  xoshiro256_state state;
  seed_xoshiro256(&state, seed_value);

  // Create R vector
  SEXP r_output = PROTECT(allocVector(REALSXP, n_values));

  // Get a pointer to the double data of the R vector
  double* vec_data = REAL(r_output);

  // Generate a random number and store it in the array
  for(int i = 0; i < n_values; i++) {
    vec_data[i] = xoshiro_uniform(&state);
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
    xoshiro256_state state;
    seed_xoshiro256(&state, seed_value);

    // Create R vector
    SEXP r_output = PROTECT(allocVector(REALSXP, n_values));

    // Get a pointer to the double data of the R vector
    double* vec_data = REAL(r_output);

    // Generate a random number and store it in the array
    for(int i = 0; i < n_values; i++) {
        vec_data[i] = (double) next(&state);
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
    xoshiro256_state state;
    seed_xoshiro256(&state, seed_value);

    // Protect the input SEXP
    PROTECT(r_vector);

    // Get a pointer to the integer data of the R vector
    int* vec_data = INTEGER(r_vector);

    // Shuffle the array using the Fisher-Yates (or Knuth shuffle) algorithm
    for (int i = vector_length - 1; i > 0; i--) {
        int j = next(&state) % (i + 1); // generates random index between 0 and i
        int tmp = vec_data[j];
        vec_data[j] = vec_data[i];
        vec_data[i] = tmp;
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
    xoshiro256_state state;
    seed_xoshiro256(&state, seed_value);

    // Create R vector
    SEXP r_output = PROTECT(allocVector(REALSXP, vector_length));

    // Get a pointer to the double data of the R vector
    double* vec_data = REAL(r_output);

    // Shuffle
    for(int i = 0; i < vector_length; i++) {
        vec_data[i] = (double) (next(&state) % vector_length) + 1;
    }

    // Release protected SEXP objects
    UNPROTECT(1);

    // Return the result
    return r_output;

}

/* References

 // xoshiro256++
 Blackman, D., & Vigna, S. (2021). Scrambled linear pseudorandom
 number generators. ACM Transactions on Mathematical Software
 (TOMS), 47(4), 1-32. https://doi.org/10.1145/3460772

 // splitmix64
 Steele Jr, G. L., Lea, D., & Flood, C. H. (2014). Fast splittable
 pseudorandom number generators. ACM SIGPLAN Notices, 49(10), 453-472.
 https://doi.org/10.1145/2714064.2660195

 // upper bound of overlapping subsequences
 Vigna, S. (2020). On the probability of overlap of random subsequences
 of pseudorandom number generators. Information Processing Letters,
 158, 105939. https://doi.org/10.1016/j.ipl.2020.105939


*/
