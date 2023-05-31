#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include "modularity.h"

SEXP r_signed_modularity(SEXP r_input_network, SEXP r_input_memberships) {

    // Obtain columns
    int cols = ncols(r_input_network);

    // Call the C functions
    struct modularity_result Q_values = modularity_values(REAL(r_input_network), cols);
    double Q = signed_modularity(Q_values, INTEGER(r_input_memberships), cols);

    // Free memory
    free(Q_values.positive_modularity_values);
    free(Q_values.negative_modularity_values);

    // Wrap the result into a SEXP object
    SEXP r_modularity = PROTECT(allocVector(REALSXP, 1));
    REAL(r_modularity)[0] = Q;
    UNPROTECT(1);

    // Return the result
    return r_modularity;

}
