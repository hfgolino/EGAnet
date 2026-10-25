/*  Written in 2019 by David Blackman and Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide.

Permission to use, copy, modify, and/or distribute this software for any
purpose with or without fee is hereby granted.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR
IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE. */
 
 // Modified by Alexander P. Christensen, 15 August 2023.

#include <stdint.h>
#include <R.h>
#include <Rinternals.h>
#include "nanotime.h"
#include "xoshiro.h"

// Process ID (used to make clock-based seeds unique across forked workers)
#ifdef _WIN32
#include <process.h>
#define get_process_id _getpid
#else
#include <unistd.h>
#define get_process_id getpid
#endif

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
  source file because they are not used in {L0ggm}

 `jump` ensures non-overlapping subsequences within xoshiro256++;
 however, {L0ggm} does not take advantage of it

 instead, seeds are pre-generated and allow splitmix64 to
 do some of the work

 according to Blackman and Vigna (2019), splitmix64 passes BigCrush,
 meaning there is sufficient randomness when setting xoshiro256++'s state

 further, in a footnote (p.11), they state that "with 256 bits of state,
 2^64 sequences of length 2^64 starting at 2^64 random points in the
 state space have an overlap probability less than 2^-64"

 leaving the sequence length aside, there is never a case that people
 should ever get close to 2^64 random starting points (seeds) in {L0ggm}
 (e.g., 500 is the default for `bootEGA`) -- even with Monte Carlo
 simulation work most values wouldn't cross 2^20 or 1 million
 random starting points

 based on Vigna (2020), the upper bound on the overlap probability
 is formally defined as n^2 * L / P where n = processors (or starting points),
 L is sequence length, and P is the period of the PRNG

 using a sequence length of 2^64 (far beyond anything used in {L0ggm}
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

// Number of clock-based seeds drawn by this process
static uint64_t seed_counter = 0;

/*

 Seed for calls that do not supply one (seed = 0)

 the clock alone is not unique: parallel workers (e.g., `mclapply`)
 that call at the same moment read the same time in nanoseconds
 and would then share an identical random stream

 the clock, process ID, and a per-process counter are therefore
 chained through splitmix64 (a bijection at each step) so that
 different workers (process ID) and different calls within a worker
 (counter) cannot share a seed, even when the clock values are equal

 R is single-threaded, so the counter does not need to be atomic

*/
static uint64_t random_seed(void) {

    // Clock time in nanoseconds
    uint64_t x = get_time_ns();
    uint64_t seed = splitmix64(&x);

    // Mix in process ID
    x = seed ^ (uint64_t) get_process_id();
    seed = splitmix64(&x);

    // Mix in call counter
    x = seed ^ ++seed_counter;
    return splitmix64(&x);

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
  if(seed_value == 0) { // Use clock time, process ID, and call counter
    seed_value = random_seed();
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
    if(seed_value == 0) { // Use clock time, process ID, and call counter
      seed_value = random_seed();
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
    if(seed_value == 0) { // Use clock time, process ID, and call counter
      seed_value = random_seed();
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

// Function to weighted shuffle indices without replacement
SEXP r_xoshiro_weighted_shuffle(SEXP r_vector, SEXP r_prob, SEXP r_seed) {

    // Get length of R vector
    int vector_length = length(r_vector);

    // Initialize seed
    uint64_t seed_value = (uint64_t) REAL(r_seed)[0];

    // For random seed, use zero
    if(seed_value == 0) { // Use clock time, process ID, and call counter
        seed_value = random_seed();
    }

    // Seed the random number generator
    xoshiro256_state state;
    seed_xoshiro256(&state, seed_value);

    // Protect inputs
    PROTECT(r_vector);
    PROTECT(r_prob);

    // Get pointers to vector and probability data
    int    *vec_data = INTEGER(r_vector);
    double *prob     = REAL(r_prob);

    // Copy probs into working array and compute total weight
    double *w = (double*) R_alloc(vector_length, sizeof(double));
    double  total = 0.0;
    for(int i = 0; i < vector_length; i++) {
        w[i]  = prob[i];
        total += prob[i];
    }

    // Weighted Fisher-Yates shuffle over full vector
    for(int i = 0; i < vector_length - 1; i++) {

        // Draw uniform on remaining total weight
        double u      = ((next(&state) >> 11) * (1.0 / (UINT64_C(1) << 53))) * total;
        double cumsum = 0.0;
        int    j      = i; // fallback to current position

        // Find the index whose cumulative weight first exceeds u
        for(int k = i; k < vector_length; k++) {
            cumsum += w[k];
            if(u <= cumsum) {
                j = k;
                break;
            }
        }

        // Swap vector elements
        int tmp_v   = vec_data[j];
        vec_data[j] = vec_data[i];
        vec_data[i] = tmp_v;

        // Swap weight elements to keep indices and weights aligned
        double tmp_w = w[j];
        w[j]   = w[i];
        w[i]   = tmp_w;

        // Remove selected weight from total
        total -= tmp_w;

    }

    // Unprotect inputs
    UNPROTECT(2);

    // Return vector
    return r_vector;

}

// Function to shuffle indices with replacement
SEXP r_xoshiro_shuffle_replace(SEXP r_vector, SEXP r_size, SEXP r_seed) {

  // Get length of R vector
  int vector_length = length(r_vector);

  // Get output size (default to vector_length if size <= 0)
  int size = asInteger(r_size);
  if (size <= 0) {
    size = vector_length;
  }

  // Initialize R values
  uint64_t seed_value = (uint64_t) REAL(r_seed)[0];

  // For random seed, use zero
  if(seed_value == 0) { // Use clock time, process ID, and call counter
    seed_value = random_seed();
  }

  // Seed the random number generator
  xoshiro256_state state;
  seed_xoshiro256(&state, seed_value);

  // Create R vector of length `size`
  SEXP r_output = PROTECT(allocVector(REALSXP, size));

  // Get a pointer to the double data of the R vector
  double* vec_data = REAL(r_output);

  // Sample with replacement
  for(int i = 0; i < size; i++) {
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
