#ifndef XOSHIRO128PLUSPLUS_H
#define XOSHIRO128PLUSPLUS_H

#include <stdint.h>

// Function to get the next random number
uint32_t next(void);

// Function to set single seed to get the 4 random seeds
void seed_xoshiro256(uint32_t seed);

// Function to generate random uniform data between 0 and 1
float xoshiro_uniform(void);

#endif /* XOSHIRO128PLUSPLUS_H */
