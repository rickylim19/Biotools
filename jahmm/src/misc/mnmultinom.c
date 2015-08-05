#include "mnmultinom.h"

int
is_invalid
(
    const int *yz,
    int k,
    int r
)
// SYNOPSIS:                                                              
//   Helper function for `mnmultinom_prob`.                               
//                                                                        
// NAs of type 'int' is the largest negative value. More generally, any   
// negative value in 'yz' is invalid.                                     
{
   for (int i = 0 ; i < r ; i++) if (yz[i + k*r] < 0) return 1;
   return 0;
}

void
mnmultinom_prob
(
   // input //
   const int *n_states,
   const int *n_obs,
   const int *dim_yz,
   const int * restrict yz,
   // params //
   const double * restrict t,
   const double * restrict a,
   const double * restrict p,
   const double * restrict q,
   // index //
         int * restrict index,
   // control //
   const int * restrict output,
   // output //
   double * restrict pem
)
// SYNOPSIS:                                                             
//   Compute emission probabilities with a mixture negative multinomial  
//   model. Since those are up to a multiplicative constant in the       
//   forward-backward algorithm, we can drop the multiplicative terms    
//   that do not depend on the state of the HMM.                         
//   Since the negative nultinomial takes discrete values, we can cache  
//   the results for reuse in order to save computation. This is done    
//   by indexing the series.                                             
//                                                                       
//   My parametrization is of the form:                                  
//                                                                       
//   theta * p_0(i)^a * p_1(i)^y * p_2(i)^z_1 * ... * p_r+1(i)^z_r +     
//   (1-theta) * q_0(i)^a * q_1(i)^y * q_2(i)^z_1 * ... * q_r+1(i)^z_r   
//                                                                       
// NUMERICAL STABILITY:                                                  
//   Each term of the sum above is computed in log space, the result is  
//   the computed as the sum of two exponentials. NA emissions are       
//   allowed and yield NA for the whole line of emissions.               
//                                                                       
// ARGUMENTS:                                                            
//   'n_states': (1) number of states in the HMM (alias 'm')             
//   'n_obs': (1) length of the sequence of observations (alias 'n')     
//   'dim_yz': (1) number of columns of 'yz' (alias 'r')                 
//   'yz': (n_obs,dim_yz) profiles                                       
//   'a': (1) model parameter                                            
//   'p': (dim_yz,m) model parameters                                    
//   'q': (dim_yz,m) model parameters                                    
//   'output': the type of output to produce (see below)                 
//   'pem': (n_obs,n_states) emission probability                        
//                                                                       
// RETURN:                                                               
//   'void'                                                              
//                                                                       
// SIDE EFFECTS:                                                         
//   Update 'pem' in place.                                              
//                                                                       
// OUTPUT:                                                               
//   The output type for 'pem' can be the emission probability in        
//   linear space (1), the same emission probability in log space (2),   
//   the ratio of probabilities of the mixture model (3), or in linear   
//   by default and in log space in case of underflow (0). 'output'      
//   also controls the verbosity. If the third bit is set, i.e. the      
//   value is set to 4, 5 or 6, the function will suppress warnings.     
//   Setting the third bit of 'output' (ie 4, 5, 6 or 7) turns off the   
//   verbosity. Setting the fourth bit of 'output' forces to compute     
//   the constant leading terms in emission probabilities.               
{

   char *depends   = "compute in lin space, log space if underflow",
        *log_space = "always compute in log space",
        *ratio     = "compute the probability ratio of the subtypes",
        *lin_space = "always compute in linear space";
   char *cases[4] = {depends, log_space, ratio, lin_space};
   char *output_type = cases[*output & 3];

   int compute_constant_terms = (*output >> 3) & 1;

   int n = *n_obs;
   int m = *n_states;
   int r = *dim_yz;

   if (*index < 0) indexts(n, r, yz, index);

   double logp[(r+1)*m];
   double logq[(r+1)*m];

   // If the third bit of 'output' is set, suppress warnings
   // by setting 'warned' to 1.
   int warned = (*output >> 2) & 1;

   // Make sure that 'p' and 'q' define probabilities.
   // If not, renormalize them.
   for (int i = 0 ; i < m ; i++) {
      double renorm_p = 0.0;
      double renorm_q = 0.0;
      for (int j = 0 ; j < r+1 ; j++) {
         renorm_p += p[j+i*(r+1)];
         renorm_q += q[j+i*(r+1)];
      }
      int p_not_normalized = fabs(renorm_p - 1.0) > DBL_EPSILON;
      int q_not_normalized = fabs(renorm_q - 1.0) > DBL_EPSILON;
      if (!warned && (p_not_normalized || q_not_normalized)) {
         fprintf(stderr, "warning: renormalizing 'p' and/or 'q'\n");
         warned = 1;
      }
      for (int j = 0 ; j < r+1 ; j++) {
         logp[j+i*(r+1)] = log(p[j+i*(r+1)] / renorm_p);
         logq[j+i*(r+1)] = log(q[j+i*(r+1)] / renorm_q);
      }
   }

   // The following variable 'row_of_na' comes in handy to write
   // full lines of NAs in the emissions.
   double *row_of_na = malloc(m * sizeof(double));
   if (row_of_na == NULL) {
      fprintf(stderr, "memory error, cannot allocate 'row_of_na'");
      return;
   }
   for (int i = 0 ; i < m ; i++) row_of_na[i] = NAN;
   // We will also need the following terms often.
   double log_theta = log(*t);
   double log_one_minus_theta = log(1-*t);

   for (int k = 0 ; k < n ; k++) {
      // Indexing allows to compute the terms only once. If the term
      // has been computed before, copy the value and move on.
      if (index[k] < k) {
         memcpy(pem + k*m, pem + index[k]*m, m * sizeof(double));
         continue;
      }

      // This is the firt occurrence of the emission in the times
      // series. We need to compute the emission probability.

      // Test for the presence of invalid/NA emissions in the row.
      if (is_invalid(yz, k, r)) {
         memcpy(pem + k*m, row_of_na, m * sizeof(double));
         continue;
      }
      
      for (int i = 0 ; i < m ; i++) {
         double p_term = log_theta + *a * logp[0+i*(r+1)];
         double q_term = log_one_minus_theta + *a * logq[0+i*(r+1)];
         for (int j = 0 ; j < r ; j++) {
            p_term += yz[j+k*r] * logp[j+1+i*(r+1)];
            q_term += yz[j+k*r] * logq[j+1+i*(r+1)];
         }

         if (output_type == ratio) {
            // The following expression is robust to underflow, it
            // should always be between 0.0 and 1.0 included.
            pem[i+k*m] = 1 / (1+exp(q_term-p_term));
         }
         else {
            // This expression looks odd, but it is more numerically
            // stable than the naive computation.
            double small = q_term;
            double big = p_term;
            if (small > big) {
               small = p_term;
               big = q_term;
            }
            pem[i+k*m] = big + log(1+exp(small-big));
         }
      }

      if (output_type == ratio) continue;

      if (compute_constant_terms) {
         double c_term = -lgamma(*a);
         double sum = *a;
         for (int j = 0 ; j < r ; j++) {
            int term = yz[j+k*r];
            sum += term;
            c_term -= lgamma(term+1);
         }
         c_term += lgamma(sum);
         for (int i = 0 ; i < m ; i++) {
            pem[i+k*m] += c_term;
         }
      }

      if (output_type == log_space) continue;

      double sum = 0.0;
      double lin[m];
      for (int i = 0 ; i < m ; i++) sum += lin[i] = exp(pem[i+k*m]);
      if (sum > 0 || output_type == lin_space) {
         memcpy(pem+k*m, lin, m * sizeof(double));
      }

   }

   free(row_of_na);
   return;

}
