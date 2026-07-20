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

// Constants in `genz_bivariate_normal` (Genz & Ge's MVBVU)
// Gauss-Legendre points/weights on [-1, 1] (only positive half stored;
// the routine exploits symmetry), used for |rho| < 0.3, < 0.75, and >= 0.75
// respectively (giving 6-, 12-, and 20-point rules)
#define GENZ_COR_MAX 0.925
#define GENZ_HK_MIN -160.0
extern const double GENZ_X6[3];
extern const double GENZ_W6[3];
extern const double GENZ_X12[6];
extern const double GENZ_W12[6];
extern const double GENZ_X20[10];
extern const double GENZ_W20[10];

// Shared constant (bivariate normal density, and `genz_bivariate_normal`)
#define TWO_PI 6.283185307179586
#define SQRT_TWO_PI 2.5066282746310002 // sqrt(TWO_PI), folded at compile time

// Constants in `newton_raphson` / `newton_raphson_arcsin`
// Fisher-scoring / exact-Newton estimation of rho (replaces Brent's method):
// converges in a handful of iterations given a reasonable starting value,
// so `NEWTON_MAX_ITER` is intentionally much smaller than Brent's `MAX_ITER`.
// A larger budget is kept in reserve for the (rarer) damped-step path near
// the +-1 boundary
#define NEWTON_MAX_ITER 40
#define NEWTON_TOL 1e-09 // tolerance on the (per-observation) score
#define NEWTON_PROB_MIN 1e-09 // avoid division by (near-)zero cell probabilities
#define NEWTON_MAX_STEP 0.3 // per-iteration step cap (rho scale, and theta scale in radians)
#define NEWTON_THETA_MAX 1.56 // keeps theta on its principal branch (rho up to ~0.9999)

#endif /* POLYCHORIC_H */
