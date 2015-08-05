#include <time.h>
#include <stdlib.h>
#include <pthread.h>
#include "utils.h"
#include "zinm.h"
#include "hmm.h"


#ifndef _PSO_HMM
#define _PSO_HMM

void
pso
(
   // input //
   const int *n_states,
   const int *n_obs,
   const int *dim_y,
   const int *y,
   // fixed params //
   const double *a,
   const double *pi,
   // start conditions //
         double *Q,
         double *p,
   // index //
         int *index,
   // output //
   double *loglikmax
);

#endif
