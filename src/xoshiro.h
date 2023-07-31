#ifndef XOSHIRO256PLUSPLUS_H
#define XOSHIRO256PLUSPLUS_H

#include <stdint.h>

// Function to get the next random number
uint64_t next(void);

// Get initial states
uint64_t splitmix64(uint64_t *x);

// Function to set single seed to get the 4 random seeds
void seed_xoshiro256(uint64_t seed);

// Function to generate random uniform data between 0 and 1
double xoshiro_uniform(void);

#endif /* XOSHIRO256PLUSPLUS_H */
