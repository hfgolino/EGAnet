#include <time.h>
#include <stdint.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include "shuffle.h"

// Get clock time in nanoseconds
uint64_t get_time_ns() {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}

SEXP r_time(SEXP r_scale) {

    // Create R vector
    SEXP r_output = PROTECT(allocVector(REALSXP, 1));

    // Get time
    switch(INTEGER(r_scale)[0]) {
    case 1:
        REAL(r_output)[0] = (double) get_time_ns();
        break;
    case 2:
        REAL(r_output)[0] = (double) get_time_ns() * 1e-03;
        break;
    case 3:
        REAL(r_output)[0] = (double) get_time_ns() * 1e-06;
        break;
    case 4:
        REAL(r_output)[0] = (double) get_time_ns() * 1e-09;
        break;
    case 5:
        REAL(r_output)[0] = (double) get_time_ns() * 1e-09 * 0.016667;
        break;
    case 6:
        REAL(r_output)[0] = (double) get_time_ns() * 1e-09 * 0.002778;
        break;
    default:
        REAL(r_output)[0] = (double) get_time_ns();
    }

    // Release protected SEXP objects
    UNPROTECT(1);

    // Return the result
    return r_output;

}

// Fisher-Yates (or Knuth) Shuffle
int* shuffle(int* arr, int cols, int seed) {

    // Make a copy of the input array
    int* copy = malloc(cols * sizeof(int));
    memcpy(copy, arr, cols * sizeof(int));

    // Initialize iterators
    int i, j, temp;

    // For random seed, use zero
    if(seed == 0) {
        // Use clocktime in nanoseconds
        uint64_t nano = get_time_ns();

        // Seed the random number generator
        srand(nano);
    } else {
        // Use user set seed
        srand(seed);
    }

    // Iterate through the array from the last element to the first
    for (i = cols - 1; i > 0; i--) {
        // Generate a random index between 0 and i (inclusive)
        j = rand() % (i + 1);

        // Swap the current element with the randomly selected element
        temp = copy[i];
        copy[i] = copy[j];
        copy[j] = temp;
    }

    // Return copy
    return copy;

}

SEXP r_shuffle(SEXP r_vector, SEXP r_seed) {

    // Get length of R vector
    int vector_length = length(r_vector);

    // Create R vector
    SEXP r_output = PROTECT(allocVector(INTSXP, vector_length));

    // Shuffle vector
    int* shuffled = shuffle(INTEGER(r_vector), vector_length, INTEGER(r_seed)[0]);

    // Copy the shuffled array to r_output
    memcpy(INTEGER(r_output), shuffled, vector_length * sizeof(int));

    // Free the shuffled array
    free(shuffled);

    // Release protected SEXP objects
    UNPROTECT(1);

    // Return the result
    return r_output;

}

// Shuffle with replacement
int* shuffle_replace(int* arr, int cols, int seed) {

    // Allocate an array for the result
    int* result = malloc(cols * sizeof(int));

    // For random seed, use zero
    if(seed == 0) {
        // Use clocktime in nanoseconds
        uint64_t nano = get_time_ns();

        // Seed the random number generator
        srand(nano);
    } else {
        // Use user set seed
        srand(seed);
    }

    // For each element in the result
    for (int i = 0; i < cols; i++) {
        // Generate a random index between 0 and cols-1 (inclusive)
        int j = rand() % cols;

        // Set the current element to a randomly selected element from the original array
        result[i] = arr[j];
    }

    // Return result
    return result;

}

SEXP r_shuffle_replace(SEXP r_vector, SEXP r_seed) {

    // Get length of R vector
    int vector_length = length(r_vector);

    // Create R vector
    SEXP r_output = PROTECT(allocVector(INTSXP, vector_length));

    // Shuffle vector
    int* shuffled_replaced = shuffle_replace(INTEGER(r_vector), vector_length, INTEGER(r_seed)[0]);

    // Copy the shuffled array to r_output
    memcpy(INTEGER(r_output), shuffled_replaced, vector_length * sizeof(int));

    // Free the shuffled array
    free(shuffled_replaced);

    // Release protected SEXP objects
    UNPROTECT(1);

    // Return the result
    return r_output;

}
