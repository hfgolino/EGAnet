#ifndef XOSHIRO256PLUSPLUS_H
#define XOSHIRO256PLUSPLUS_H

#include <stdint.h>

// Structure for seeding states
typedef struct {
    uint64_t s[4];
} xoshiro256_state;


// Function to get the next random number
uint64_t next(xoshiro256_state* state);

// Get initial states
uint64_t splitmix64(uint64_t *x);

// Function to set single seed to get the 4 random seeds
void seed_xoshiro256(xoshiro256_state* state, uint64_t seed);

// Function to generate random uniform data between 0 and 1
double xoshiro_uniform(xoshiro256_state* state);

#endif /* XOSHIRO256PLUSPLUS_H */
