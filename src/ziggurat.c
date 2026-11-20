// Original code by John Burkardt, distributed under the GNU LGPL license.
// Modified by Alexander P. Christensen, 30 July 2023.

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "ziggurat.h"
#include "xoshiro.h"
#include "nanotime.h"

// Seed from R (xoshiro.c)
uint64_t seed_from_r(SEXP r_seed);

// Set up static values for normal table
static uint32_t kn[128];
static double fn[128];
static double wn[128];
static int initialized = 0;

// Function to initialize table
void r4_nor_initialize(void) {
  if (!initialized) {
    r4_nor_setup(kn, fn, wn);
    initialized = 1;
  }
}

/*

 The base ziggurat.c code has been modified to support the use
 of the xoshiro256++ generator for random uniform numbers. The
 xoshiro256++ uniform generator replaced the `r4_uni`,
 which was the original random uniform generating function.
 This new uniform generating function is called `xoshiro_uniform`.

 This replacement voids the need for several functions that
 were previously used to set random seeds:
 +   `cong_seeded`
 +   `kiss_seeded`
 +   `mwc_seeded`
 +   `shr3_seeded`

 This replacement eliminates the need to call a seed into `r4_uni` as
 well as any manipulation of the seed. Code for the xoshiro256++
 can be found in the xoshiro.c file.

 Further, all `float` have been updated to `double` to be compatible
 with the 64-bit generation used in xshiro256++

 Modified date: 30.07.2023
 Modified by: Alexander P. Christensen <alexpaulchristensen@gmail.com>

 */

/******************************************************************************/

double r4_nor ( xoshiro256_state* state, uint32_t kn[128], double fn[128], double wn[128] )

  /******************************************************************************/
  /*
   Purpose:

   r4_nor() returns a normally distributed single precision real value.

   Discussion:

   The value returned is generated from a distribution with mean 0 and
   variance 1.

   The underlying algorithm is the ziggurat method.

   Before the first call to this function, the user must call R4_NOR_SETUP
   to determine the values of KN, FN and WN.

   Thanks to Chad Wagner, 21 July 2014, for noticing a bug of the form
   if ( x * x <= y * y );   <-- Stray semicolon!
   {
   break;
   }

   Licensing:

   This code is distributed under the GNU LGPL license.

   Modified:

   21 July 2014

   Author:

   John Burkardt

   Reference:

   George Marsaglia, Wai Wan Tsang,
   The Ziggurat Method for Generating Random Variables,
   Journal of Statistical Software,
   Volume 5, Number 8, October 2000, seven pages.

   Parameters:

   Input/output, uint32_t *JSR, the seed.

   Input, uint32_t KN[128], data computed by R4_NOR_SETUP.

   Input, float FN[128], WN[128], data computed by R4_NOR_SETUP.

   Output, float R4_NOR, a normally distributed random value.
   */
{
  /*
   Modified 20.11.2026: each step draws 64 bits; the signed 32-bit value
   `hz` is taken from the upper 32 bits and the strip index `iz` from the
   lowest 7 bits, so that the index and the value are independent (the
   original used the same low bits for both; see Doornik, 2005). Absolute
   values are computed in 64 bits (no overflow for the most negative value)
  */
  int32_t hz;
  uint32_t iz;
  uint64_t bits;
  double x, y, value;

  bits = next ( state );
  hz = ( int32_t ) ( uint32_t ) ( bits >> 32 );
  iz = ( uint32_t ) ( bits & 127 );

  if ( llabs ( ( int64_t ) hz ) < ( int64_t ) kn[iz] ) {
    value = ( double ) ( hz ) * wn[iz];
  } else {

    for ( ; ; ) {

      if ( iz == 0 ) {

        for ( ; ; ) {

          x = - 0.2904764 * log ( xoshiro_uniform( state ) );
          y = - log ( xoshiro_uniform( state ) );
          if ( x * x <= y + y ) {
            break;
          }

        }

        value = ( hz <= 0 ) ? ( -ZIGGURAT_R - x ) : ( ZIGGURAT_R + x );

        break;

      }

      x = ( double ) ( hz ) * wn[iz];

      if (
          fn[iz] + xoshiro_uniform( state ) *
          ( fn[iz-1] - fn[iz] ) <
          exp ( - 0.5 * x * x )
      ) {
        value = x;
        break;
      }

      bits = next ( state );
      hz = ( int32_t ) ( uint32_t ) ( bits >> 32 );
      iz = ( uint32_t ) ( bits & 127 );

      if ( llabs ( ( int64_t ) hz ) < ( int64_t ) kn[iz] ) {
        value = ( double ) ( hz ) * wn[iz];
        break;
      }

    }

  }

  return value;

}
/******************************************************************************/

void r4_nor_setup ( uint32_t kn[128], double fn[128], double wn[128] )

  /******************************************************************************/
  /*
   Purpose:

   r4_nor_setup() sets data needed by r4_nor().

   Licensing:

   This code is distributed under the GNU LGPL license.

   Modified:

   14 October 2013

   Author:

   John Burkardt

   Reference:

   George Marsaglia, Wai Wan Tsang,
   The Ziggurat Method for Generating Random Variables,
   Journal of Statistical Software,
   Volume 5, Number 8, October 2000, seven pages.

   Parameters:

   Output, uint32_t KN[128], data needed by R4_NOR.

   Output, float FN[128], WN[128], data needed by R4_NOR.
   */
{
  
  /*
    
   Original code for comparison is provide below.
   Constants were calculated by hand computing each
   calculation
    
   
      double dn = 3.442619855899;
      int i;
      const double m1 = 2147483648.0;
      double q;
      double tn = 3.442619855899;
      const double vn = 9.91256303526217E-03;
    
      q = vn / exp ( - 0.5 * dn * dn );
    
      kn[0] = ( uint32_t ) ( ( dn / q ) * m1 );
      kn[1] = 0;
    
      wn[0] = ( q / m1 );
      wn[127] = ( dn / m1 );
    
      fn[0] = 1.0;
      fn[127] = ( exp ( - 0.5 * dn * dn ) );
    
      for ( i = 126; 1 <= i; i-- )
      {
        dn = sqrt ( - 2.0 * log ( vn / dn + exp ( - 0.5 * dn * dn ) ) );
        kn[i+1] = ( uint32_t ) ( ( dn / tn ) * m1 );
        tn = dn;
        fn[i] = ( exp ( - 0.5 * dn * dn ) );
        wn[i] = ( dn / m1 );
      }
  */
  
  int i;
  double dn = DN;
  double tn = DN;

  kn[0] = 1991057938;
  kn[1] = 0;

  wn[0] = 1.729040521542796563654E-09;
  wn[127] = 1.603094793809112161398E-09;

  fn[0] = 1.0;
  fn[127] = 0.002669629083880925288704;

  for ( i = 126; 1 <= i; i-- )
  {
    dn = sqrt ( - 2.0 * log ( VN / dn + exp ( - 0.5 * dn * dn ) ) );
    kn[i+1] = ( uint32_t ) ( ( dn / tn ) * M1 );
    tn = dn;
    fn[i] = ( exp ( - 0.5 * dn * dn ) );
    wn[i] = ( dn / M1 );
  }

  return;
}

/******************************************************************************/

SEXP r_ziggurat(SEXP n, SEXP r_seed) {

  // Initialize values
  int n_values = INTEGER(n)[0];

  // Seed (0 = clock time, process ID, and call counter, as in xoshiro.c)
  uint64_t seed_value = seed_from_r(r_seed);

  // Seed the (uniform) random number generator
  xoshiro256_state state;
  seed_xoshiro256(&state, seed_value);

  // Initialize table (if necessary)
  r4_nor_initialize();

  // Create R vector
  SEXP r_output = PROTECT(allocVector(REALSXP, n_values));
  
  // Get a pointer to the double data of the R vector
  double* vec_data = REAL(r_output);
  
  // Generate random numbers
  for(int i = 0; i < n_values; i++) {
    vec_data[i] = r4_nor(&state, kn, fn, wn);
  }

  // Release protected SEXP objects
  UNPROTECT(1);

  // Return result
  return r_output;

}
