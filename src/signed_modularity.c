#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include "modularity.h"

SEXP r_signed_modularity(SEXP r_input_network, SEXP r_input_memberships, SEXP r_resolution) {

    // Obtain columns
    int cols = ncols(r_input_network);

    // Call the C functions
    struct modularity_result Q_values = modularity_values(
      REAL(r_input_network), cols, REAL(r_resolution)[0]
    );

    // Initialize R object
    SEXP r_modularity = PROTECT(allocVector(REALSXP, 1));

    // Calculate signed modularity
    REAL(r_modularity)[0] = signed_modularity(Q_values, INTEGER(r_input_memberships), cols);

    // Free R output
    UNPROTECT(1);

    // Free memory
    free(Q_values.positive_modularity_values);
    free(Q_values.negative_modularity_values);

    // Return the result
    return r_modularity;

}
