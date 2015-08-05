#include <emmintrin.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

#ifndef _PUBLIC_ZINM_HEADER
#define _PUBLIC_ZINM_HEADER

#define HISTO_INIT_SIZE 128
#define ZINM_MAXITER 32
#define ZINM_TOL 1e-6


#define sq(x) ((x)*(x))

struct zinb_part_t;

typedef struct zinb_par_t zinb_par_t;

struct zinb_par_t {
   double   a;
   double   pi;
   double   p;
};

double       eval_nb_dfda(double, const tab_t *);
double       eval_nb_f(double, const tab_t *);
double 	    eval_zinb_dfda(double, double, unsigned int);
double 	    eval_zinb_dfdp(double, double, unsigned int, double);
double 	    eval_zinb_dgda(double, double, const tab_t *);
double 	    eval_zinb_f(double, double, unsigned int, double);
double 	    eval_zinb_g(double, double, const tab_t *);
double	    ll_zinb(double, double, double, const tab_t *);
zinb_par_t * mle_nb(int *, unsigned int);
zinb_par_t * mle_zinb(int *, unsigned int);
double       nb_est_alpha(tab_t *);
zinb_par_t * new_zinb_par(size_t);
#endif
