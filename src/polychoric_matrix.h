#ifndef POLYCHORIC_H
#define POLYCHORIC_H

// Constant for hard cut-off for polychoric
#define CUT 11 // similar to {Turbofuns}

// Constants in `bsm_inverse_cdf`
extern const double CONST_A[6];
extern const double CONST_B[5];
extern const double CONST_C[6];
extern const double CONST_D[4];

// Constants in `joint_frequency_table`
#define MISSING 99

// Constants in `error_function`
#define A1 0.254829592
#define A2 -0.284496736
#define A3 1.421413741
#define A4 -1.453152027
#define A5 1.061405429
#define P 0.3275911

// Constants in `drezner_bivariate_normal`
#define INT_NX 5
#define COR_MAX 0.7
#define BV_FAC1 0.13298076
#define BV_FAC2 0.053051647
extern const double DOUBLE_X[5];
extern const double DOUBLE_W[5];

// Constants in `optimize`
#define LOWER -1.0
#define UPPER 1.0
#define TOL 1e-05 // tolerate to floating point
#define MAX_ITER 100
#define ZEPS 1e-10

#endif /* POLYCHORIC_H */
