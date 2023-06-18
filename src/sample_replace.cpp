#include <Rcpp.h>
#include <random>

// Sample with replacement
int sample_with_replacement(int n, std::mt19937& g) {
    std::uniform_int_distribution<> dis(1, n);
    return dis(g); // generates a random index between 1 and n
}

//' @noRd
// [[Rcpp::export]]
Rcpp::IntegerVector r_sample_with_replacement(int n, unsigned int seed) {
    // Initialize random number generator
    std::mt19937 g(seed);

    // Initialize vector for R
    Rcpp::IntegerVector r_with_replacement(n);

    // Generate numbers
    for (int i = 0; i < n; i++) {
        r_with_replacement[i] = sample_with_replacement(n, g);
    }

    return r_with_replacement;
}
