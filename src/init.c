#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Declare the C functions you want to make available to R here
extern SEXP r_signed_louvain(SEXP r_input_network);
extern SEXP r_signed_modularity(SEXP r_input_network, SEXP r_input_memberships);
extern SEXP r_polychoric_correlation_matrix(SEXP r_input_matrix, SEXP r_empty_method, SEXP r_empty_value);
extern SEXP _EGAnet_r_sample_with_replacement(SEXP n, SEXP seed);
extern SEXP _EGAnet_r_sample_without_replacement(SEXP arr, SEXP seed);
extern SEXP _EGAnet_r_sample_seeds(SEXP n, SEXP seed);
extern SEXP _EGAnet_r_random_uniform(SEXP n, SEXP min, SEXP max, SEXP seed);


// Register native routine
static const R_CallMethodDef CallEntries[] = {

    {
        "r_signed_louvain", // Name of function call in R
        (DL_FUNC)&r_signed_louvain, // Name of C function
         1 // Number of arguments
    },
    {
        "r_signed_modularity", // Name of function call in R
        (DL_FUNC)&r_signed_modularity, // Name of C function
         2 // Number of arguments
    },
    {
        "r_polychoric_correlation_matrix", // Name of function call in R
        (DL_FUNC)&r_polychoric_correlation_matrix, // Name of C function
         3 // Number of arguments
    },
    {
        "_EGAnet_r_sample_with_replacement", // Name of function call in R
        (DL_FUNC)&_EGAnet_r_sample_with_replacement, // Name of C function
         2 // Number of arguments
    },
    {
        "_EGAnet_r_sample_without_replacement", // Name of function call in R
        (DL_FUNC)&_EGAnet_r_sample_without_replacement, // Name of C function
         2 // Number of arguments
    },
    {
        "_EGAnet_r_sample_seeds", // Name of function call in R
        (DL_FUNC)&_EGAnet_r_sample_seeds, // Name of C function
         2 // Number of arguments
    },
    {
        "_EGAnet_r_random_uniform", // Name of function call in R
        (DL_FUNC)&_EGAnet_r_random_uniform, // Name of C function
         4 // Number of arguments
    },
    {NULL, NULL, 0}

};

// Set up call in package
void R_init_EGAnet(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
