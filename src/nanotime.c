#include <time.h>
#include <stdint.h>
#include <R.h>
#include <Rinternals.h>

// Get clock time in nanoseconds
uint64_t get_time_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}
