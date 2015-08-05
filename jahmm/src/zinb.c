#include "zinb.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


double
eval_nb_f
(
         double a,
   const tab_t *tab
)
{

   double retval;
   double prev;
   // Convenience variables.
   const unsigned int *val = tab->val;
   const unsigned int *num = tab->num;

   size_t nobs = num[0];
   double mean = num[0] * val[0];

   prev = digamma(a + val[0]);
   retval = num[0] * prev;
   // Iterate over the occurrences and compute the new value
   // of digamma either by the recurrence relation, or by
   // a new call to 'digamma()', whichever is faster.
   for (size_t i = 1 ; i < tab->size ; i++) {
      nobs += num[i];
      mean += num[i] * val[i];
      prev = (val[i] - val[i-1] == 1) ?
         prev + 1.0 / (a-1 + val[i]) :
         digamma(a + val[i]);
      retval += num[i] * prev;
   }

   mean /= nobs;
   retval += nobs*(log(a) - digamma(a) - log(a + mean));

   return retval;

}


double
eval_nb_dfda
(
   double a,
   const tab_t *tab
)
{

   double retval;
   double prev;
   // Convenience variables.
   const unsigned int *val = tab->val;
   const unsigned int *num = tab->num;

   size_t nobs = num[0];
   double mean = num[0] * val[0];

   prev = trigamma(a + val[0]);
   retval = num[0] * prev;
   // Iterate over the occurrences and compute the new value
   // of trigamma either by the recurrence relation, or by
   // a new call to 'trigamma()', whichever is faster.
   for (size_t i = 1 ; i < tab->size ; i++) {
      nobs += num[i];
      mean += num[i] * val[i];
      prev = (val[i] - val[i-1] == 1) ?
         prev - 1.0 / sq(a-1 +val[i]) :
         trigamma(a + val[i]);
      retval += num[i] * prev;
   }

   mean /= nobs;
   retval += nobs*(mean/(a*(a+mean)) - trigamma(a));

   return retval;

}


double
eval_zinb_f
(
   double a,
   double p,
   unsigned int   nz,
   double sum
)
{
   return a*nz / (p*(1-pow(p,a))) - sum/(1-p);
}


double
eval_zinb_g
(
         double   a,
         double   p,
   const tab_t  * tab
)
{

   // Convenience variables.
   const unsigned int *val = tab->val;
   const unsigned int *num = tab->num;

   unsigned int nz = 0;
   double retval = 0.0;
   double prev = digamma(a + val[0]);

   if (val[0] > 0) {
      retval += num[0] * prev;
      nz += num[0];
   }

   // Iterate over the occurrences and compute the new value
   // of digamma either by the recurrence relation, or by
   // a new call to 'digamma()', whichever is faster.
   const size_t imin = val[0] == 0 ? 1 : 0;
   for (size_t i = imin ; i < tab->size ; i++) {
      nz += num[i];
      prev = (val[i] - val[i-1] == 1) ?
         prev + 1.0 / (a-1 + val[i]) :
         digamma(a + val[i]);
      retval += num[i] * prev;
   }

   retval += nz*(log(p) / (1-pow(p,a)) - digamma(a));
   return retval;

}


double
eval_zinb_dfda
(
   double a,
   double p,
   unsigned int   nz
)
{
   const double ppa = pow(p,a);
   return nz*(1-ppa+a*ppa*log(p)) / (p*sq(1-ppa));
}

double
eval_zinb_dfdp
(
   double a,
   double p,
   unsigned int   nz,
   double m
)
{
   const double ppa = pow(p,a);
   return -(nz*a*(1-(a+1)*ppa) / sq(p*(1-ppa)) + m/sq(1-p));
}


double
eval_zinb_dgda
(
         double   a,
         double   p,
   const tab_t  * tab
)
{
   
   // Convenience variables.
   const unsigned int *val = tab->val;
   const unsigned int *num = tab->num;
   const double ppa = pow(p,a);

   unsigned int nz = 0;
   double retval = 0.0;
   double prev = trigamma(a + val[0]);

   if (val[0] > 0) {
      retval += num[0] * prev;
      nz += num[0];
   }

   // Iterate over the occurrences and compute the new value
   // of digamma either by the recurrence relation, or by
   // a new call to 'trigamma()', whichever is faster.
   const size_t imin = val[0] == 0 ? 1 : 0;
   for (size_t i = imin ; i < tab->size ; i++) {
      nz += num[i];
      prev = (val[i] - val[i-1] == 1) ?
         prev - 1.0 / sq(a-1 + val[i]) :
         trigamma(a + val[i]);
      retval += num[i] * prev;
   }

   retval += nz*(sq(log(p))*ppa / sq(1-ppa) - trigamma(a));
   return retval;

}


double
ll_zinb
(
         double   a,
         double   p,
         double   pi,
   const tab_t  * tab
)
{

   // Convenience variables.
   const unsigned int *val = tab->val;
   const unsigned int *num = tab->num;
   const unsigned int z0 = val[0] == 0 ? num[0] : 0;
   const double logp_ = log(1-p);

   unsigned int nobs = z0;
   double retval = z0*log(pi*pow(p,a) + 1-pi);

   const size_t imin = val[0] == 0 ? 1 : 0;
   for (size_t i = imin ; i < tab->size ; i++) {
      retval += num[i] * (lgamma(a+val[i]) + val[i]*logp_);
      nobs += num[i];
   }
   // Ignore constant factorial terms.
   retval += (nobs-z0) * (a*log(p) - lgamma(a) + log(pi));

   return retval;

}


zinb_par_t *
mle_zinb
(
   int *x,
   unsigned int nobs
)
{

   tab_t *tab = tabulate(x, nobs);

   double sum = 0.0;
   unsigned int nona = 0;
   for (size_t i = 0 ; i < tab->size ; i++) {
      sum += tab->val[i]*tab->num[i];
      nona += tab->num[i];
   }

   // Extract the number of all-zero observaions.
   const unsigned int z0 = tab->val[0] == 0 ? tab->num[0] : 0;

   double deficit[11] = {0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1};
   double init_a[12] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1};
   double init_p[12] = {0,0,0,0,0,0,0,0,0,0,0,.5};

   // Deplete some 0s from the observations, compute alpha
   // and p0 with standard negative binomial estimation
   // and keep the values to be used as initial conditions.
   for (size_t i = 0 ; i < 11 ; i++) {
      if (tab->val[0] == 0) tab->num[0] = z0 * (1-deficit[i]);
      double newmean = sum / nona / (1.0 - z0*deficit[i]/nobs);
      double alpha = nb_est_alpha(tab);
      init_a[i] = alpha;
      init_p[i] = alpha / (alpha + newmean);
   }

   // Reset 'tab'.
   if (tab->val[0] == 0) tab->num[0] = z0;

   zinb_par_t *par = calloc(1, sizeof(zinb_par_t));
   if (par == NULL) {
      fprintf(stderr, "memory error: %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }

   // Try initial conditions. Number 12 is a safety in case
   // all the rest failed during the first phase.
   double max_loglik = -1.0/0.0;
   for (size_t i = 0 ; i < 12 ; i++) {

      if (init_a[i] < 0) continue;  // Skip failures.

      double a = init_a[i];
      double p = init_p[i];

      double grad;
      unsigned int iter = 0;

      double f = eval_zinb_f(a, p, nobs-z0, sum);
      double g = eval_zinb_g(a, p, tab);

      // Newton-Raphson iterations.
      while ((grad = f*f+g*g) > sq(ZINM_TOL) && iter++ < ZINM_MAXITER) {

         double dfda, dfdp, dgda, dgdp;
         dfda = dgdp = eval_zinb_dfda(a, p, nobs-z0);
         dfdp = eval_zinb_dfdp(a, p, nobs-z0, sum);
         dgda = eval_zinb_dgda(a, p, tab);

         double denom = dfdp*dgda - dfda*dgdp;
         double da = (f*dgdp - g*dfdp) / denom;
         double dp = (g*dfda - f*dgda) / denom;
         // Maintain 'a' and 'p' in their domain of definition.
         while (a+da < 0 || p+dp < 0 || p+dp > 1) {
            da /= 2;
            dp /= 2;
         }
         f = eval_zinb_f(a+da, p+dp, nobs-z0, sum);
         g = eval_zinb_g(a+da, p+dp, tab);
         // Backtrack if necessary.
         for (int j = 0 ; j < ZINM_MAXITER && f*f+g*g > grad ; j++) {
            da /= 2;
            dp /= 2;
            f = eval_zinb_f(a+da, p+dp, nobs-z0, sum);
            g = eval_zinb_g(a+da, p+dp, tab);
         }

         a = a+da;
         p = p+dp;

      }

      double pi = (nobs-z0) / (1-pow(p,a)) / nobs;
      if (pi > 1) pi = 1.0;
      if (pi < 0) pi = 0.0;
      double loglik = ll_zinb(a, p, pi, tab);
      if (loglik > max_loglik) {
         max_loglik = loglik;
         par->a = a;
         par->pi = pi;
         par->p = p;
      }
            
   }

   free(tab);

   return par;

}


double
nb_est_alpha
(
   tab_t *tab
)
{

   // Find upper and lower bouds for a(lpha).
   double a = 1.0;
   double a_lo;
   double a_hi;
   if (eval_nb_f(a, tab) < 0) {
      a /= 2;
      while (eval_nb_f(a, tab) < 0) a /= 2;
      a_lo = a;
      a_hi = a*2;
   }
   else {
      a *= 2;
      while (eval_nb_f(a, tab) > 0) a *= 2;
      a_lo = a/2;
      a_hi = a;
   }

   // Input is pathological.
   if (a_lo > 128) return -1.0;

   double new_a = (a_lo + a_hi) / 2;
   for (int i = 0 ; i < ZINM_MAXITER ; i++) {
      a = (new_a < a_lo || new_a > a_hi) ?
         (a_lo + a_hi) / 2 :
         new_a;
      double f = eval_nb_f(a, tab);
      if (f < 0) a_hi = a; else a_lo = a;
      if ((a_hi - a_lo) < ZINM_TOL) break;
      double dfda = eval_nb_dfda(a, tab);
      new_a = a - f / dfda;
   }

   return a;

}

zinb_par_t *
mle_nb
(
   int *x,
   unsigned int nobs
)
{

   tab_t *tab = tabulate(x, nobs);
   double a = nb_est_alpha(tab);

   // Return NULL if failed.
   if (a < 0) return NULL;

   // Allocate return struct.
   zinb_par_t *par = calloc(1, sizeof(zinb_par_t));
   if (par == NULL) {
      fprintf(stderr, "memory error: %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }

   double sum = 0.0;
   unsigned int nona = 0;
   for (size_t i = 0 ; i < tab->size ; i++) {
      nona += tab->num[i];
      sum += tab->val[i]*tab->num[i];
   }

   free(tab);

   par->a = a;
   par->p = a / (a + (sum/nona));

   return par;

}
