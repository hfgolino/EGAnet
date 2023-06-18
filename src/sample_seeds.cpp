#include <Rcpp.h>
#include <random>

//' @noRd
// [[Rcpp::export]]
Rcpp::IntegerVector r_sample_seeds(int n, unsigned int seed = 0) {
    
    // Initialize random number generator
    std::mt19937 g(seed == 0 ? std::random_device{}() : seed);

    // Uniform distribution
    std::uniform_int_distribution<> dist(0, 2147483647); // 32-bit machine maximum

    // Generate numbers
    Rcpp::IntegerVector r_seeds(n);
    for (int i = 0; i < n; i++) {
        r_seeds[i] = dist(g);
    }

    return r_seeds;
}
