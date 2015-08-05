#include "jahmm.h"

void
apologize
(void)
{
   char msg[] =
   "  We apologize for the trouble you are experiencing.\n"
   "  Please contact guillaume.filion@gmail.com about\n"
   "  this issue and attach the error message above.\n";
   fprintf(stderr, "%s", msg);
}


unsigned int
nobs
(
   const ChIP_t *ChIP
)
{
   unsigned int retval = 0;
   for (size_t i = 0 ; i < ChIP->nb ; i++) retval += ChIP->size[i];
   return retval;
}


jahmm_t *
do_jahmm
(
   ChIP_t *ChIP
)
{

   // Extract the dimensions of the observations.
   const unsigned int r = ChIP->r;
   const unsigned int n = nobs(ChIP);


   // Extract the first ChIP profile, which is the sum of
   // negative controls.
   int *ctrl = malloc(n * sizeof(int));
   if (ctrl == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }

   for (size_t i = 0 ; i < n ; i++) {
      ctrl[i] = ChIP->y[0+i*r];
   }

   zinb_par_t *z = mle_zinb(ctrl, n);
   if (z == NULL) {
      fprintf(stderr, "jahmm failure %s:%d\n", __FILE__, __LINE__);
      apologize();
      return NULL;
   }

   free(ctrl);

   // Jahmm uses 3 states.
   const unsigned int m = 3;

   double *Q = malloc(m*m * sizeof(double));
   double *p = malloc(m*(r+1) * sizeof(double));
   if (Q == NULL || p == NULL) {
      fprintf(stderr, "memory error: %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }
   // Set initial values of 'Q'.
   for (size_t i = 0 ; i < m ; i++) {
   for (size_t j = 0 ; j < m ; j++) {
      Q[i+j*m] = (i==j) ? .95 : .05/(m-1);
   }
   }

   // Set initial values of 'p'. The are not normalize, but the
   // call to 'bw_zinm' will normalize them.
   for (size_t i = 0 ; i < m ; i++) {
      p[0+i*(r+1)] = z->p;
      p[1+i*(r+1)] = 1 - z->p;
      for (size_t j = 2 ; j < r+1 ; j++) {
         p[j+i*(r+1)] = i + 0.5;
      }
   }

   jahmm_t *jahmm = new_jahmm(m, ChIP);
   if (jahmm == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      free(Q);
      free(p);
      return NULL;
   }
   set_jahmm_par(jahmm, Q, z->a, z->pi, p);

   free(Q);
   free(p);

   Q = jahmm->Q;
   p = jahmm->p;

   // Run the Baum-Welch algorithm.
   bw_zinm(jahmm);

   // Reorder the states in case they got scrambled.
   // We use the value of p0 as a sort key.
   int tmp, map[3] = {0,1,2};
   if (p[0*(r+1)] < p[1*(r+1)]) { tmp = map[0]; map[0] = map[1]; map[1] = tmp; }
   if (p[1*(r+1)] < p[2*(r+1)]) { tmp = map[1]; map[1] = map[2]; map[2] = tmp; }
   if (p[0*(r+1)] < p[1*(r+1)]) { tmp = map[0]; map[0] = map[1]; map[1] = tmp; }

   if (map[0] != 0 || map[1] != 1) {
      // The states have beem scrambled. We need to reorder
      // 'Q', 'p', 'phi' and 'pem'.
      
      double *phi = jahmm->phi;
      double *pem = jahmm->pem;

      double *Q_ = malloc(m*m * sizeof(double));
      double *p_ = malloc(m*(r+1) * sizeof(double));
      double *phi_ = malloc(m*n * sizeof(double));
      double *pem_ = malloc(m*n * sizeof(double));
      if (Q_ == NULL || p_ == NULL || phi_ == NULL || pem_ == NULL) {
         // TODO: free everything
         fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
         return NULL;
      }

      for (size_t j = 0 ; j < m ; j++) {
      for (size_t i = 0 ; i < m ; i++) {
         Q_[map[i]+map[j]*m] = Q[i+j*m];
      }
      }
      memcpy(jahmm->Q, Q_, m*m * sizeof(double));

      for (size_t j = 0 ; j < m ; j++) {
      for (size_t i = 0 ; i < r+1 ; i++) {
         p_[i+map[j]*(r+1)] = p[i+j*(r+1)];
      }
      }
      memcpy(jahmm->p, p_, m*(r+1) * sizeof(double));

      for (size_t j = 0 ; j < m ; j++) {
      for (size_t i = 0 ; i < n ; i++) {
         phi_[map[j]+i*m] = phi[j+i*m];
         pem_[map[j]+i*m] = pem[j+i*m];
      }
      }
      memcpy(jahmm->phi, phi_, m*n * sizeof(double));
      memcpy(jahmm->pem, pem_, m*n * sizeof(double));

      free(Q_);
      free(p_);
      free(phi_);
      free(pem_);

   }

   // Run the Viterbi algorithm.
   int *path = malloc(n * sizeof(int));
   double *initp = malloc(m * sizeof(double));
   double *log_Q = malloc(m*m * sizeof(double));
   if (path == NULL || initp == NULL || log_Q == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      // TODO: free everything.
      return NULL;
   }

   for (size_t i = 0 ; i < m ; i++) initp[i] = log(1.0/m);
   for (size_t i = 0 ; i < m*m ; i++) log_Q[i] = log(Q[i]);

   block_viterbi(m, ChIP->nb, ChIP->size, log_Q, initp, jahmm->pem, path);
   jahmm->path = path;

   free(initp);
   free(log_Q);

   return jahmm;

}


int
is_invalid
(
    const int * y,
          int   k,
          int   r
)
// SYNOPSIS:                                                              
//   Helper function. NAs of type 'int' is the largest negative
//   value. More generally, any negative value in 'y' is invalid.
{
   for (int i = 0 ; i < r ; i++) if (y[i + k*r] < 0) return 1;
   return 0;
}

int
is_all_zero
(
   const int * y,
         int   k,
         int   r
)
// SYNOPSIS:                                                              
//   Helper function for `zinm_prob`. Returns 1 if and only if all
//   the observations are 0.
{
   for (int i = 0 ; i < r ; i++) if (y[i + k*r] != 0) return 0;
   return 1;
}


void
zinm_prob
(
         jahmm_t * restrict jahmm,
   const int     * restrict index,
   // control //
         int                otype,
   // output //
         double  * restrict pem
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
//        p_0(i)^a * p_1(i)^y_1 * p_2(i)^y_2 * ... * p_r+1(i)^y_r        
//                                                                       
//   And in the case that all emissions are 0                            
//                                                                       
//                      pi * p_0(i)^a + (1-pi)                           
//                                                                       
// NUMERICAL STABILITY:                                                  
//   Each term of the sum above is computed in log space, the result is  
//   the computed as the sum of two exponentials. NA emissions are       
//   allowed and yield NA for the whole line of emissions.               
//                                                                       
// ARGUMENTS:                                                            
//   'ChIP': struct of observations                                      
//   'par': struct of parameters for the ZINM distribution               
//   'index': a precomputed index of the ChIP data                       
//   'otype': the type of output to produce (see below)                  
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
//   log space (1), the same emission probability in linear space (2),   
//   or in linear by default and in log space in case of underflow (0).  
//   'otype' also controls the verbosity. If the third bit is set,       
//   i.e. the value is set to 4, 5 or 6, the function will suppress      
//   warnings. Setting the fourth bit of 'otype' forces to compute       
//   the constant terms in emission probabilities.                       
{

   ChIP_t *ChIP = jahmm->ChIP;
   unsigned int temp = 0;
   for (size_t i = 0 ; i < ChIP->nb ; i++) {
      temp += ChIP->size[i];
   }

   const unsigned int   r  = ChIP->r;
   const int          * y  = ChIP->y;
   const unsigned int   m  = jahmm->m;
   const double         a  = jahmm->a;
   const double         pi = jahmm->pi;
   const double       * p  = jahmm->p;
   const unsigned int   n  = temp;

   char *depends   = "compute in lin space, log space if underflow";
   char *log_space = "always compute in log space";
   char *lin_space = "always compute in linear space";
   char *cases[3] = {depends, log_space, lin_space};
   char *output_type = cases[otype & 3];

   int compute_constant_terms = (otype >> 3) & 1;

   double *logp = malloc((r+1)*m * sizeof(double));
   if (logp == NULL) {
      fprintf(stderr, "memory error (%s:%d)\n", __FILE__, __LINE__);
      return;
   }

   // If the third bit of 'otype' is set, suppress warnings
   // by setting 'warned' to 1.
   int warned = (otype >> 2) & 1;

   // Make sure that 'p' defines a probability.
   for (size_t i = 0 ; i < m ; i++) {
      double sump = 0.0;
      for (size_t j = 0 ; j < r+1 ; j++) {
         // Cannot normalize negative values. Sorry folks.
         if (p[j+i*(r+1)] < 0) {
            fprintf(stderr, "error: 'p' negative\n");
            return;
         }
         sump += p[j+i*(r+1)];
      }
      int p_normalized_no = fabs(sump - 1.0) > DBL_EPSILON;
      if (!warned && p_normalized_no) {
         fprintf(stderr, "warning: renormalizing 'p'\n");
         warned = 1;
      }
      for (int j = 0 ; j < r+1 ; j++) {
         logp[j+i*(r+1)] = log(p[j+i*(r+1)] / sump);
      }
   }

   // The following variable 'row_of_na' comes in handy to write
   // full lines of NAs in the emissions.
   double *row_of_na = malloc(m * sizeof(double));
   if (row_of_na == NULL) {
      fprintf(stderr, "memory error (%s:%d)\n", __FILE__, __LINE__);
      return;
   }
   for (int i = 0 ; i < m ; i++) row_of_na[i] = NAN;

   for (int k = 0 ; k < n ; k++) {
      // Indexing allows to compute the terms only once. If the term
      // has been computed before, copy the value and move on.
      if (index[k] < k) {
         // TODO: bypass the cache for writing. This would save
         // some time when tere are a lot of observations.
         // This would require reformattnig 'pem' and using
         // _mm_stream_pd(double *p, __m128d a)
         memcpy(pem + k*m, pem + index[k]*m, m * sizeof(double));
         continue;
      }

      // This is the firt occurrence of the emission in the times
      // series. We need to compute the emission probability.
      // Test the presence of invalid/NA emissions in the row.
      // If so, fill the row with NAs and move on.
      if (is_invalid(y, k, r)) {
         memcpy(pem + k*m, row_of_na, m * sizeof(double));
         continue;
      }

      if (is_all_zero(y, k, r)) {
         // Emissions are all zeros, use the zero-inflated
         // term from the zinm model.
         for (int i = 0 ; i < m ; i++) {
            pem[i+k*m] = log(pi*exp(a*logp[0+i*(r+1)]) + (1.0-pi));
         }
      }
      else {
         // Otherwise use the standard probability.
         for (int i = 0 ; i < m ; i++) {
            pem[i+k*m] = a * logp[0+i*(r+1)];
            for (int j = 0 ; j < r ; j++) {
               pem[i+k*m] += y[j+k*r] * logp[(j+1)+i*(r+1)];
            }
         }
      }

      if (compute_constant_terms) {
         double c_term = -lgamma(a);
         double sum = a;
         for (int j = 0 ; j < r ; j++) {
            int term = y[j+k*r];
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

   free(logp);
   free(row_of_na);
   return;

}


ChIP_t *
read_file
(
   FILE *inputf
)
{

   ssize_t nread;
   size_t nchar = 256; 
   char *line = malloc(256 * sizeof(char));
   if (line == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }

   // Count the lines of the file.
   size_t nobs = 0;
   while ((nread = getline(&line, &nchar, inputf)) != -1) nobs++;
   rewind(inputf);
   nobs--; // Discount the header.

   // Read and parse the header separately.
   nread = getline(&line, &nchar, inputf);
   if (nread == -1) return NULL;

   // Use the header to get the number of tokens 'ntok'.

   char *s;
   int k;

   const char sep = '\t';
   for ( k = 0, s = line ; (s = strchr(s, sep)) != NULL ; s++) k++;

   const size_t dim = k;

   int *y = malloc(nobs*dim * sizeof(int));
   if (y == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }

   histo_t *histo = new_histo();
   size_t current = -1;
   char prevblk[256] = {0};
   unsigned int lineno = 0;

   for ( ; (nread = getline(&line, &nchar, inputf)) != -1 ; lineno++) {

      // Trim new line character.
      if (line[nread-1] == '\n') line[nread-1] = '\0';

      // Update block if needed.
      s = strchr(line, sep);
      if (s == NULL) {
         fprintf(stderr, "error (1) parsing line %d:\n%s\n",
               lineno+2, line);
         return NULL;
      }
      *s = '\0';
      if (strcmp(line, prevblk) != 0) {
         strncpy(prevblk, line, 255);
         current++;
      }

      histo_push(&histo, current);

      // Gather values.
      char *e;
      for ( k = 0, s++ ; (e = strchr(s, sep)) != NULL ; k++, s = e+1) {
         if (k > dim) {
            fprintf(stderr, "too many entries line %d\n", lineno+2);
            return NULL;
         }
         *e = '\0';
         if (strcmp(s, "NA") == 0) {
            y[k + lineno*dim] = -1;
         }
         else {
            int v = (int) strtol(s, &s, 10);
            if (s != e) {
               fprintf(stderr, "non integer entry line %d\n", lineno+2);
               return NULL;
            }
            y[k + lineno*dim] = v;
         }
      }

      // Last token.
      if (strcmp(s, "NA") == 0) {
         y[k + lineno*dim] = -1;
      }
      else {
         int v = (int) strtol(s, &s, 10);
         if (*s != '\0') {
            fprintf(stderr, "non integer entry line %d\n", lineno+2);
            return NULL;
         }
         y[k + lineno*dim] = v;
      }
      if (k != (dim-1)) {
         fprintf(stderr, "too few entries line %d\n", lineno+2);
         return NULL;
      }

   }

   free(line);

   tab_t *tab = compress_histo(histo);
   ChIP_t *ChIP = new_ChIP(dim, tab->size, y, tab->num);

   free(histo);
   free(tab);

   if (ChIP == NULL) return NULL;
   return ChIP;

}


void
update_trans
(
         size_t   m,
         double * Q,
   const double * trans
)
{

   for (size_t i = 0 ; i < m ; i++) {
      double sum = 0.0;
      for (size_t j = 0 ; j < m ; j++) {
         sum += trans[i+j*m];
      }
      for (size_t j = 0 ; j < m ; j++) {
         Q[i+j*m] = trans[i+j*m] / sum;
      }
   }

   return;

}


double
eval_bw_f
(
   double a,
   double pi,
   double p0,
   double A,
   double B,
   double C,
   double D,
   double E
)
{
   double term1 = (D + a*A) / p0;
   double term2 = B * pi*a*pow(p0,a-1) / (pi*pow(p0,a)+1-pi);
   return p0 + E/(term1 + term2) - 1.0 / C;
}

double
eval_bw_dfdp0
(
   double a,
   double pi,
   double p0,
   double A,
   double B,
   double C,
   double D,
   double E
)
{

   double term1 = (D + a*A) / p0;
   double term2 = B * pi*a*pow(p0,a-1) / (pi*pow(p0,a)+1-pi);
   double subterm3a = (1-pi)*pi*a*(a-1)*pow(p0,a-2);
   double subterm3b = sq(pi)*a*pow(p0,2*a-2);
   double term3 = B * (subterm3a - subterm3b) / sq(pi*pow(p0,a)+1-pi);
   double term4 = (D + a*A) / sq(p0);

   return 1 - E/sq(term1 + term2) * (term3-term4);

}


void
bw_zinm
(
   jahmm_t *jahmm   
)
{

   // Unpack parameters.
   ChIP_t *ChIP = jahmm->ChIP;
   size_t temp = 0;
   for (size_t i = 0 ; i < ChIP->nb ; i++) {
      temp += ChIP->size[i];
   }

   // Constants.
   const size_t         n    = temp;
   const size_t         m    = jahmm->m;
   const size_t         r    = ChIP->r;
   const unsigned int   nb   = ChIP->nb;
   const unsigned int * size = ChIP->size;
   const int          * y    = ChIP->y;
   const double         a    = jahmm->a;
   const double         pi   = jahmm->pi;
   const double         R    = (jahmm->p[1]) / jahmm->p[0];

   // Variables optimized by the Baum-Welch algorithm.
   double *p = jahmm->p;
   double *Q = jahmm->Q;

   // Check the input.
   for (size_t i = 1 ; i < m ; i++) {
      double ratio = jahmm->p[1+i*(r+1)] /  jahmm->p[0+i*(r+1)];
      if (fabs(ratio - R) > 1e-3) {
         fprintf(stderr, "warning (%s): 'p' inconsistent\n", __func__);
      }
   }

   int *index = malloc(n * sizeof(int));
   double *pem = malloc(n*m * sizeof(double));
   double *phi = malloc(n*m * sizeof(double));
   if (index == NULL || pem == NULL || phi == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      return;
   }

   double *trans = malloc(m*m * sizeof(double));
   double *ystar = malloc(r * sizeof(double));
   double *prob = malloc(m * sizeof(double));
   if (trans == NULL || ystar == NULL || prob == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      return;
   }

   double *newp = malloc(m*(r+1) * sizeof(double));

   for (size_t i = 0 ; i < m ; i++) prob[i] = 1.0 / m;

   // Index the time series now. This would be done by
   // 'zinm_prob()' anyway, but we will need the index of
   // the first all-0 emission later.
   int i0 = indexts(n, r, y, index);

   // Start Baum-Welch cycles.
   for (jahmm->iter = 1 ; jahmm->iter < BW_MAXITER ; jahmm->iter++) {

#ifndef NDEBUG
fprintf(stderr, "iter: %d\r", jahmm->iter);
#endif

      // Update emission probabilities and run the block
      // forward backward algorithm.
      unsigned int lin_space_no_warn = 4;
      zinm_prob(jahmm, index, lin_space_no_warn, pem);
      jahmm->l = block_fwdb(m, nb, size, Q, prob, pem, phi, trans);

      // Update 'Q'.
      update_trans(m, Q, trans);

      // Update 'p'.
      for (size_t i = 0 ; i < m ; i++) {
         // Compute the constants.
         double A = 0.0;
         double B = 0.0;
         double C = 1+R;
         double D = 0.0;
         double E = 0.0;
         memset(ystar, 0, r * sizeof(double));
         for (size_t k = 0 ; k < n; k++) {
            if (index[k] == i0) {
               B += phi[i+k*m];
            }
            else {
               // Skip invalid entries.
               if (is_invalid(y, k, r)) continue;
               A += phi[i+k*m];
               D += phi[i+k*m] * y[0+k*r];
               for (size_t j = 1 ; j < r ; j++) {
                  ystar[j] += phi[i+k*m] * y[j+k*r];
               }
            }
         }
         for (size_t j = 1 ; j < r ; j++) {
            E += ystar[j];
         }
         
         // Find upper and lower bound for 'p0'.
         double p0 = .5;
         double p0_lo;
         double p0_hi;
         if (eval_bw_f(a, pi, p0, A, B, C, D, E) < 0) {
            p0 *= 2;
            while (eval_bw_f(a, pi, p0, A, B, C, D, E) < 0) p0 *= 2;
            p0_lo = p0/2;
            p0_hi = p0;
         }   
         else {
            p0 /= 2;
            while (eval_bw_f(a, pi, p0, A, B, C, D, E) > 0) p0 /= 2;
            p0_lo = p0;
            p0_hi = p0*2;
         }   

         if (p0_lo > 1.0 || p0_hi < 0.0) {
            fprintf(stderr, "cannot complete Baum-Welch algorithm\n");
            free(index);
            free(newp);
            free(pem);
            free(phi);
            free(prob);
            free(trans);
            free(ystar);
            return;
         }

         double new_p0 = (p0_lo + p0_hi) / 2;
         for (int j = 0 ; j < BT_MAXITER ; j++) {
            p0 = (new_p0 < p0_lo || new_p0 > p0_hi) ?
               (p0_lo + p0_hi) / 2 : 
               new_p0;
            double f = eval_bw_f(a, pi, p0, A, B, C, D, E);
            if (f > 0) p0_hi = p0; else p0_lo = p0;
            if ((p0_hi - p0_lo) < TOLERANCE) break;
            double dfdp0 = eval_bw_dfdp0(a, pi, p0, A, B, C, D, E);
            new_p0 = p0 - f / dfdp0;
         }   

         // Update the state-independent parameters.
         newp[0+i*(r+1)] = p0;
         newp[1+i*(r+1)] = p0 * R;
         // Now udpaate the state-dependent parameters.
         double term1 = (D + a*A) / p0;
         double term2 = B * pi*a*pow(p0,a-1) / (pi*pow(p0,a)+1-pi);
         double normconst = (term1 + term2) / C;
         for (size_t j = 1 ; j < r ; j++) {
            newp[(j+1)+i*(r+1)] = ystar[j] / normconst;
         }

      }
      
      // Check convergence
      double maxd = 0.0;
      for (size_t i = 0 ; i < m*(r+1) ; i++) {
         double thisd = fabs(newp[i]-p[i]);
         maxd = thisd > maxd ? thisd : maxd;
      }

      if (maxd < TOLERANCE) break;
      memcpy(p, newp, m*(r+1) * sizeof(double));

   }

#ifndef NDEBUG
fprintf(stderr, "\n");
#endif

   free(newp);
   free(prob);
   free(trans);
   free(ystar);

   // Compute final emission probs in log space.
   unsigned int log_space_no_warn = 5;
   zinm_prob(jahmm, index, log_space_no_warn, pem);

   free(index);

   // 'Q','p' and 'l' have been updated in-place.
   jahmm->phi = phi;
   jahmm->pem = pem;

   return;

}


jahmm_t *
new_jahmm
(
   unsigned int   m,
   ChIP_t       * ChIP
)
{

   if (ChIP == NULL) return NULL;
   unsigned int r = ChIP->r;

   jahmm_t *new = calloc(1, sizeof(jahmm_t));
   double *newQ = malloc(m*m * sizeof(double));
   double *newp = malloc(m*(r+1) * sizeof(double));
   if (new == NULL || newQ == NULL || newp == NULL) return NULL;

   new->m = m;
   new->ChIP = ChIP;
   new->Q = newQ;
   new->p = newp;

   return new;

}


void
set_jahmm_par
(
   jahmm_t      * jahmm,
   double const * Q, 
   double         a,
   double         pi,
   double const * p
)
{

   const unsigned int m = jahmm->m;
   const unsigned int r = jahmm->ChIP->r;
   memcpy(jahmm->Q, Q, m*m * sizeof(double));
   memcpy(jahmm->p, p, m*(r+1) * sizeof(double));
   jahmm->a = a;
   jahmm->pi = pi;

   return;

}


ChIP_t *
new_ChIP
(
         unsigned int   r,
         unsigned int   nb,
                  int * y,
   const unsigned int * size
)
{

   size_t extra = nb * sizeof(unsigned int);
   ChIP_t *new = calloc(1, sizeof(ChIP_t) + extra);
   if (new == NULL) return NULL;

   new->r = r;
   new->nb = nb;
   new->y = y;
   memcpy(new->size, size, nb * sizeof(unsigned int));

   return new;

}

void
destroy_jahmm_all
(
   jahmm_t * jahmm
)
{

   if (jahmm->ChIP != NULL) {
      if (jahmm->ChIP->y != NULL) free(jahmm->ChIP->y);
      free(jahmm->ChIP);
   }
   if (jahmm->Q != NULL) free(jahmm->Q);
   if (jahmm->p != NULL) free(jahmm->p);
   if (jahmm->phi != NULL) free(jahmm->phi);
   if (jahmm->pem != NULL) free(jahmm->pem);
   if (jahmm->path != NULL) free(jahmm->path);
   free(jahmm);

   return;

}
