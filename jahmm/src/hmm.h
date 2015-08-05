#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <err.h>

#ifndef _HMM_HEADER_
#define _HMM_HEADER_

#define D double
#define I int
#define U unsigned int
#define V void
#define cD const double
#define cU const unsigned int

// function      ( 1    2   3    4    5    6   7   8  )
D  block_fwdb    (  U,  U, cU*,  D*,  D*,  D*, D*, D* );
I  block_viterbi ( cU, cU, cU*, cD*, cD*, cD*, I*     );
V  bwd           (  U,  U, cD*,  D*,  D*,  D*         );
D  fwd           (  U,  U, cD*, cD*,  D*              );
D  fwdb          (  U,  U, cD*, cD*,  D*,  D*, D*     );
V  viterbi       (  U,  U, cD*, cD*,  cD*, I*         );

#undef D
#undef I
#undef U
#undef V
#undef cD
#undef cU

#endif
