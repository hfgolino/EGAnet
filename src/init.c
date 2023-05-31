#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Declare the C functions you want to make available to R here
extern SEXP r_signed_louvain(SEXP r_input_network);
extern SEXP r_signed_modularity(SEXP r_input_network, SEXP r_input_memberships);

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
    {NULL, NULL, 0}

};

// Set up call in package
void R_init_EGAnet(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
