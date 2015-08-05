#define _GNU_SOURCE
#include <stdlib.h>
#include <time.h>
#include <glib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include "utils.h"
#include "mnmultinom.h"
#include "hmm.h"
#include "pso.h"

// --           extra declarations needed for testing           -- //

// FROM: pso.c
typedef struct {
   double l;      // log-likelihood
   double a;      // alpha
   double t;      // theta
   double *Q;
   double *p;
   double *q;
   int m;
   int r;
} params;

params *
params_change
(
   params *new,
   params *old,
   int which
);

params *
params_new
(
   int m,
   int r
);

params *
params_set
(
   params *par,
   double t,
   double a,
   double *Q,
   double *p,
   double *q
);

params *
params_cpy
(
   params *dest,
   params *par
);

void params_destroy (params *p);

double
loglik
(
   const int n,
         params *par,
   const int *yz, 
         int *index,
         double *pem 
);

