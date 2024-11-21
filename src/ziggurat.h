#include "xoshiro.h"

// Constants
#define R 3.442620 // r
#define DN 3.442619855899 // dn
#define M1 2147483648.0 // m1
#define VN 9.91256303526217E-03 // vn

// Function prototypes
void r4_nor_setup ( uint32_t kn[128], double fn[128], double wn[128] );
double r4_nor ( xoshiro256_state* state, uint32_t kn[128], double fn[128], double wn[128] );
