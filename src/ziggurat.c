#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "ziggurat.h"
#include "xoshiro.h"
#include "nanotime.h"

// Set up static values for normal table
static uint32_t kn[128];
static double fn[128];
static double wn[128];
static int initialized = 0;

// Function to initialize table
void r4_nor_initialize() {
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

double r4_nor ( uint32_t kn[128], double fn[128], double wn[128] )
  
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
  int hz;
  uint32_t iz;
  const double r = 3.442620;
  double value;
  double x;
  double y;
  
  hz = ( int ) next ( ); // cast may cause negative (and that's OK)
  iz = ( hz & 127 );
  
  if ( abs ( hz ) < kn[iz] )
  {
    value = ( double ) ( hz ) * wn[iz];
  }
  else
  {
    for ( ; ; )
    {
      if ( iz == 0 )
      {
        for ( ; ; )
        {
          x = - 0.2904764 * log ( xoshiro_uniform( ) );
          y = - log ( xoshiro_uniform( ) );
          if ( x * x <= y + y )
          {
            break;
          }
        }
        
        if ( hz <= 0 )
        {
          value = - r - x;
        }
        else
        {
          value = + r + x;
        }
        break;
      }
      
      x = ( double ) ( hz ) * wn[iz];
      
      if ( fn[iz] + xoshiro_uniform( ) * ( fn[iz-1] - fn[iz] )
             < exp ( - 0.5 * x * x ) )
      {
        value = x;
        break;
      }
      
      hz = ( int ) next ( ); // cast may cause negative (and that's OK)
      iz = ( hz & 127 );
      
      if ( abs ( hz ) < kn[iz] )
      {
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
  
  return;
}

/******************************************************************************/

SEXP r_ziggurat(SEXP n, SEXP r_seed) {
  
  // Initialize iterators
  int i;
  
  // Initialize values
  int n_values = INTEGER(n)[0];
  uint64_t seed_value = (uint64_t) REAL(r_seed)[0];
  
  // For random seed, use zero
  if(seed_value == 0) {
    
    // Use clocktime in nanoseconds
    seed_value = get_time_ns();
    
  }
  
  // Seed the (uniform) random number generator
  seed_xoshiro256(seed_value);
  
  // Initialize table (if necessary)
  r4_nor_initialize();
  
  // Create R vector
  SEXP r_output = PROTECT(allocVector(REALSXP, n_values));
  
  // Generate random numbers
  for(i = 0; i < n_values; i++) {
    REAL(r_output)[i] = r4_nor(kn, fn, wn);
  }
  
  // Release protected SEXP objects
  UNPROTECT(1);
  
  // Return result
  return r_output;
  
}
