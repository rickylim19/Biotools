#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "utils.h"

#ifndef _MIXTURE_NEG_BINOM
#define _MIXTURE_NEG_BINOM
void
mnmultinom_prob
(
   // input //
   const int *n_states,
   const int *n_obs,
   const int *dim_yz,
   const int *yz,
   // params //
   const double *t,
   const double *a,
   const double *p,
   const double *q,
   // index //
         int *index,
   // control //
   const int *output,
   // output //
   double *pem
);
#endif
