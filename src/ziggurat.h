#include "xoshiro.h"

void r4_nor_setup ( uint32_t kn[128], double fn[128], double wn[128] );
double r4_nor ( xoshiro256_state* state, uint32_t kn[128], double fn[128], double wn[128] );
