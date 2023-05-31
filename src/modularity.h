#ifndef MODULARITY_H
#define MODULARITY_H

#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>

// Structure for `modularity_values`
struct modularity_result {
    double* positive_modularity_values;
    double* negative_modularity_values;
    double positive_sum_flag;
    double negative_sum_flag;
    double positive_total_sum;
    double negative_total_sum;
};

// Function prototypes
struct modularity_result modularity_values(double* network, int cols);
double signed_modularity(struct modularity_result Q_values, int* membership, int cols);

#endif /* MODULARITY_H */
