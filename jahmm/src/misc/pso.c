#include "pso.h"

#define ITER 250
#define RAND (1+(0.5*(drand48()-.5)))

// Global lock for mutex.
pthread_mutex_t lock;

// Basic parameters.
typedef struct {
   double l;      // log-likelihood
   double a;      // alpha
   double pi;     // pi
   double *Q;
   double *p;
   int m;
   int r;
} params_t;

// Arguments of `particle`.
typedef struct {
         int n;
         params_t *opt;
   const int *y;
         int *index;
         int *signal;
} parg;


// --                  memory functions                  -- //

params_t *
params_new
(
   int m,
   int r
)
{
   params_t *par = malloc(sizeof(params_t));
   if (par == NULL) return NULL;
   // Allocate a single array for 'Q', 'p' and 'q'.
   double *array = malloc((m*(m+(r+1))) * sizeof(double));
   if (array == NULL) {
      free(par);
      return NULL;
   }
   par->Q = array;
   par->p = array + m*m;
   par->m = m;
   par->r = r;

   return par;

}


params_t *
params_set
(
   params_t *par,
   double pi,
   double a,
   double * restrict Q,
   double * restrict p
)
{

   int m = par->m;
   int r = par->r;
   par->a = a;
   par->pi = pi;
   memcpy(par->Q, Q, m*m * sizeof(double));
   memcpy(par->p, p, m*(r+1) * sizeof(double));

   return par;

}


params_t *
params_cpy
(
   params_t *dest,
   params_t *par
)
{

   int m = dest->m;
   int r = dest->r;
   dest->l = par->l;
   dest->a = par->a;
   dest->pi = par->pi;
   memcpy(dest->Q, par->Q, m*(m+(r+1)) * sizeof(double));

   return dest;

}


void
params_destroy
(params_t *p)
{
   // 'Q', 'p' and 'q' are allocated together and 'Q'
   // points to the beginning of the array.
   free(p->Q);
   free(p);
}


// --                  main functions                  -- //

double
loglik
(
   const int n,
         params_t * restrict par,
   const int * restrict y,
         int * restrict index,
         // output //
         double * restrict pem
)
{

   int m = par->m;
   int r = par->r;
   // Turn off verbosity of `mnmultinom_prob`. Otherwise irrelevant
   // messages about the normalization of 'p' and 'q' may appear.
   int o = 4;
   zinm_prob(
      &m,
      &n,
      &r,
      y,
      &par->pi,
      &par->a,
      par->p,
      index,
      &o,
      pem
   );

   double init[m]; for (int i = 0 ; i < m ; i++) init[i] = 1.0/m;
   return par->l = fwd(m, n, par->Q, init, pem);

}


params_t *
params_change
(
   params_t * restrict new,
   params_t * restrict old,
   int which
)
{

   int r = new->r;
   int m = new->m;

   switch (which) {
   case 0:
      for (int i = 0 ; i < m ; i++) {
         double sumQ = 0.0;
         for (int j = 0 ; j < m ; j++) {
            sumQ += new->Q[i+j*m] = old->Q[i+j*m] * RAND;
         }
         for (int j = 0 ; j < m ; j++) new->Q[i+j*m] /= sumQ;
      }
      break;
   case 1:
      for (int i = 0 ; i < m ; i++) {
         double sump = 0.0;
         for (int j = 1 ; j < r+1 ; j++) {
            sump += new->p[j+i*(r+1)] = old->p[j+i*(r+1)] * RAND;
         }
         sump += new->p[0+i*(r+1)] = old->p[0+i*(r+1)] *
            new->p[1+i*(r+1)] / old->p[1+i*(r+1)];
         for (int j = 0 ; j < r+1 ; j++) new->p[j+i*(r+1)] /= sump;
      }
   }

   return new;

}


void *
particle
(void *arg)
{
   // Unpack arguments.
   parg *unpack = (parg *) arg;
         int n = unpack->n;
         params_t *opt = unpack->opt;
   const int *y = unpack->y;
         int *index = unpack->index;
         int *signal = unpack->signal;

   // Allocate extra set of parameters.
   double *pem = malloc(n*opt->m * sizeof(double));
   if (pem == NULL) {
      fprintf(stderr, "memory error: cannot allocate 'pem'");
      return NULL;
   }
   params_t *par = params_new(opt->m, opt->r);
   params_t *try = params_new(opt->m, opt->r);
   
   if ((par == NULL) || (try == NULL)) {
      fprintf(stderr, "memory error: cannot allocate 'par' or 'try'");
      free(pem);
      return NULL;
   }

   // Copy best conditions.
   params_cpy(par, opt);
   params_cpy(try, opt);

   // Start simulated annealing.

   double new = par->l;
   double total = 0.0;
   double T = 0.0;
   for (int iter = 0 ; iter < ITER ; iter++) {
      double old = new;
      if (*signal > 0) {
         // Copy best conditions.
         pthread_mutex_lock(&lock);
            params_cpy(par, opt);
            (*signal)--;
         pthread_mutex_unlock(&lock);
      }
      params_change(try, par, iter % 3);

      // Keep the global best hit. 
      new = loglik(n, try, y, index, pem);
      total += abs(new-old);
      if (new > opt->l) {
         if (!(new < 0)) continue;
         pthread_mutex_lock(&lock);
            *signal = 6;
            params_cpy(opt, try);
         pthread_mutex_unlock(&lock);
      }

      // Auto-adjust temperature every 25 cycles.
      if (iter % 25 == 0) {
         T = total / 25 / 0.69;
         if (T < 0.001) T = 0.001;
         total = 0.0;
      }
      else {
         T /= .95;
      }

      // Metropolis move.
      double delta = (new-old)/T;
      if ((delta > 0) || (drand48() < exp(delta))) {
         params_cpy(par, try);
      }

   }

   free(pem);
   params_destroy(par);
   params_destroy(try);

   return NULL;

}


void
pso
(
   // input //
   const int *n_states,
   const int *n_obs,
   const int *dim_y,
   const int * restrict y,
   // fixed params //
   const double * restrict pi,
   const double * restrict a,
   // start conditions //
         double * restrict Q,
         double * restrict p,
   // index //
         int * restrict index,
   // output //
   double * restrict loglikmax
)
// SYNOPSIS:
{

   //outf = fopen("tracks.txt", "w");
   int n_threads = 12;

   int m = *n_states;
   int n = *n_obs;
   int r = *dim_y;

   // Initial parameters.
   double *pem = malloc(n*m * sizeof(double));
   if (pem == NULL) {
      fprintf(stderr, "memory error: cannot allocate 'pem'");
      return;
   }
   params_t *opt = params_new(m, r);
   if (opt == NULL) {
      fprintf(stderr, "memory error: cannot allocate 'opt'");
      return;
   }

   params_set(opt, *pi, *a, Q, p);

   // Initialize random generator.
   srand48(time(NULL));
   loglik(n, opt, y, index, pem);
   free(pem);

   // Initialize mutex.
   int err = pthread_mutex_init(&lock, NULL);
   if (err) {
      fprintf(stderr, "error initializing mutex (%d)\n", err);
      return;
   }

   // Initialize arguments.
   int signal = 0;
   parg arg = {
      .n = n,
      .opt = opt,
      .y = y,
      .index = index,
      .signal = &signal,
   };

   // Allocate 'tid'.
   pthread_t *tid = (pthread_t *) malloc(n_threads * sizeof(pthread_t));
   for (int i = 0 ; i < n_threads ; i++) {
      err = pthread_create(tid+i, NULL, &particle, &arg);
      if (err) {
         fprintf(stderr, "error creating thread (%d)\n", err);
         return;
      }
   }

   // Wait for threads to return.
   for (int i = 0 ; i < n_threads ; i++) {
      pthread_join(tid[i], NULL);
   }

   memcpy(Q, opt->Q, m*m * sizeof(double));
   memcpy(p, opt->p, m*(r+1) * sizeof(double));
   *loglikmax = opt->l;

   pthread_mutex_destroy(&lock);
   params_destroy(opt);

   free(tid);

}
