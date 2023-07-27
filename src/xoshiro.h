#ifndef XOSHIRO256_H
#define XOSHIRO256_H

#include <stdint.h>

// Function to get the next random number
uint32_t next(void);

// Function to set single seed to get the 4 random seeds
void seed_xoshiro256(uint32_t seed);

#endif /* XOSHIRO256_H */
