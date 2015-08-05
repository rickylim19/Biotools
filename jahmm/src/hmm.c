#include "hmm.h"

double
fwd
(
   // input //
         unsigned int            m,
         unsigned int            n,
   const double       * restrict Q,
   const double       * restrict init,
   // output //
         double       * restrict prob
)
// SYNOPSIS:                                                             
//   Forward algorithm.                                                  
//                                                                       
// NUMERIC ROBUSTNESS:                                                   
//   This implementation is robust to NAs and to underflow. In case of   
//   NA or underflow it will ignore the emission and treat the           
//   observation as missing (i.e. only the transitions at that position  
//   will contribute to the output). If the emission probabilities are   
//   passed as negative numbers, the algorithm will assume that they are 
//   given in log space.                                                 
//                                                                       
// ARGUMENTS:                                                            
//   'm': the number of states                                           
//   'n': the length of the sequence of observations                     
//   'Q': (m,m) transition matrix ('Q[i+j*m]' is a ij transtition).      
//   'init': (m) initial probabilities                                   
//   'prob': (m,n) emission probabilities                                
//                                                                       
// RETURN:                                                               
//   The total log-likelihood.                                           
//                                                                       
// SIDE EFFECTS:                                                         
//   Replaces 'prob' by forward alphas.                                  
{

   int i;           // State index.
   int j;           // State index.
   int k;           // Position of the time series.
   double tmp[m];   // Intermediates.
   double a[m];     // Current normalized alpha.

   // Normalization constant, also used to return log-likelihood.
   double c;
   double loglik = 0.0;

   for (k = 0 ; k < n ; k++) {
      // This is an easy pattern for the branch predictor.
      if (k == 0) {
         memcpy(tmp, init, m * sizeof(double));
      }
      else {
         memset(tmp, 0, m * sizeof(double));
         for (j = 0 ; j < m ; j++) {
         for (i = 0 ; i < m ; i++) {
            tmp[j] += a[i] * Q[i+j*m];
         }
         }
      }

      // Test for missing emission probabilities.
      int na_found = 0;
      for (j = 0 ; j < m ; j++) {
         if (prob[j+k*m] != prob[j+k*m]) {
            // NA found. Ignore emissions, and update 'prob'
            // with the value of 'a'.
            memcpy(a, tmp, m * sizeof(double));
            memcpy(prob + k*m, tmp, m * sizeof(double));
            na_found = 1;
            break;
         }
      }
      // Move on if NAs were found.
      if (na_found) continue;

      c = 0.0;
      // Test if emission probabilities have underflowed.
      // NB: we use the convention that in case all 'm' emission
      // probabilies underflow, their log is returned instead. If the
      // first one is negative, they are all computed in log space.
      if (prob[0+k*m] < 0) {
         // Use an alternative computation to obviate underflow.
         // The is is slower because of the call to the function `exp`.
         // First I find the max emission probability, then I divide 'c'
         // by the exp of that value and compensate by adding the value
         // to 'loglik' directly.
         int w = 0;
         for (j = 1 ; j < m ; j++) if (prob[j+k*m] > prob[w+k*m]) w = j;
         for (j = 0 ; j < m ; j++) {
            c += a[j] = tmp[j] * exp(prob[j+k*m] - prob[w+k*m]);
         }
         // To the exception of the correction below, the rest
         // of the computation is identical.
         loglik += prob[w+k*m];
      }
      else {
         // No underflow. Continue the forward algorithm the usual way.
         for (j = 0 ; j < m ; j++) {
            c += a[j] = tmp[j] * prob[j+k*m];
         }
      }
      if (!(c > 0)) {
         // Underflow can theoretically still happen, for instance if
         // the transition to the state with highest emission probability
         // is impossible. In this (hopeless) treat the emissions as
         // missing.
         memcpy(a, tmp, m * sizeof(double));
         memcpy(prob + k*m, tmp, m * sizeof(double));
      }
      else {
         for (j = 0 ; j < m ; j++) a[j] /= c;
         memcpy(prob +k*m, a, m * sizeof(double));
         loglik += log(c);
      }

   }

   return loglik;

}


void
bwd
(
   // input //
         unsigned int            m,
         unsigned int            n,
   const double       * restrict Q,
   // output //
         double       * restrict alpha,
         double       * restrict phi,
         double       * restrict T
)
// SYNOPSIS:                                                             
//   Backward algorithm with Markovian backward smoothing.               
//                                                                       
// ARGUMENTS:                                                            
//   'm': the number of states                                           
//   'n': the length of the sequence of observations                     
//   'Q': (m,m) transition matrix ('Q[i+j*m]' is a ij transtition).      
//   'alpha': (m,n) the forward alpha probabilities                      
//   'phi': (m,n) probabilities of states given observations             
//   'T': (m,m) sum of conditional transitions probabilties              
//                                                                       
// RETURN:                                                               
//   'void'                                                              
//                                                                       
// SIDE EFFECTS:                                                         
//   Updates 'phi' and 'T' in place.                                     

{

   int     i;       // State index.
   int     j;       // State index.
   int     k;       // Position of the time series.
   double  x;       // Sum used for computation intermediates.
   double *R;       // Reverse kernel.

   R = malloc(m*m * sizeof(double));
   if (R == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      return;
   }

   // 'T[i+j*m]' is the sum of transition probabilities from state 'i'
   // to state 'j' (congruent with 'Q') conditional on the observations.
   memset(T, 0.0, m*m * sizeof(double));

   // First iteration of the backward pass.
   memset(phi, 0.0, m*n * sizeof(double));
   memcpy(phi+(n-1)*m, alpha+(n-1)*m, m * sizeof(double));

//-----------------------------------------------------------------------
// Here we work out the local reverse kernel.                            
// ak(i) is the probability of being in state i at step k given Y1,...,k 
// bk(i) is the beta function defined by the forward-backward decom-     
// position, such that ak(i)bk(i) is the probability of being in state i 
// at step k given Y1,...,n.                                             
//                                                                       
// 'R[j+i*m]' = P(Xk=i|Xk+1=j,Y1,...,n) =                                
//       P(Xk=i,Xk+1=j,Y1,...,n) / P(Xk+1=j,Y1,...,n) =                  
//       ak(i)Q(i,j)gk+1(j)bk+1(j) / ak+1(j)bk+1(j) =                    
//       ak(i)Q(i,j)gk+1(j) / sum_i ak(i)Q(i,j)gk+1(j)                   
//       ak(i)Q(i,j) / sum_i ak(i)Q(i,j)                                 
//                                                                       
// P(Xk=i|Y1,...,n) = sum_j P(Xk=i|Xk+1=j,Y1,...,n) * P(Xk+1=j|Y1,...,n) 
// which gives the line 'phi[j+k*m] += phi[i+(k+1)*m] * R[i+j*m]'.       
//-----------------------------------------------------------------------

   // Next iterations of the backward pass.
   for (k = n-2 ; k >= 0 ; k--) {
      for (j = 0 ; j < m ; j++) {
         x = 0.0;
         // Note the double assignment in the following line.
         for (i = 0 ; i < m ; i++) x += R[j+i*m] = alpha[i+k*m]*Q[i+j*m];
         for (i = 0 ; i < m ; i++) R[j+i*m] /= x;
      }
      for (j = 0 ; j < m ; j++) {
      for (i = 0 ; i < m ; i++) {
         // Use the reverse kernel to update 'phi' and 'T'.
         x = phi[i+(k+1)*m] * R[i+j*m];
         phi[j+k*m] += x;
         T[j+i*m] += x;
      }
      }
   }

   free(R);

   return;

}


double
fwdb
(
   // input //
         unsigned int            m,
         unsigned int            n,
   const double       * restrict Q,
   const double       * restrict init,
   // output //
         double       * restrict prob,
         double       * restrict phi,
         double       * restrict T
)
// SYNOPSIS:                                                             
//   Forward-backward algorithm with Markovian backward smoothing.       
//                                                                       
// NUMERIC ROBUSTNESS:                                                   
//   This implementation is robust to NAs and to underflow. In case of   
//   NA or underflow it will ignore the emission and treat the           
//   observation as missing (i.e. only the transitions at that position  
//   will contribute to the output).                                     
//                                                                       
// ARGUMENTS:                                                            
//   'm': the number of states                                           
//   'n': the length of the sequence of observations                     
//   'Q': (m,m) transition matrix ('Q[i+j*m]' is a ij transtition).      
//   'init': (m) initial probabilities                                   
//   'prob': (m,n) emission probabilities                                
//   'phi': (m,n) probabilities given observations                       
//   'T': (m,m) sum of conditional transitions probabilties              
//                                                                       
// RETURN:                                                               
//   The total log-likelihood.                                           
//                                                                       
// SIDE EFFECTS:                                                         
//   Replaces 'prob' by alphas, updates 'phi' and 'T' in place.          
{

   double loglik = fwd(m, n, Q, init, prob);
   bwd(m, n, Q, prob, phi, T);

   return loglik;
}


void
viterbi(
   // input //
         unsigned int            m,
         unsigned int            n,
   const double       * restrict log_Q,
   const double       * restrict log_i,
   const double       * restrict log_p,
   // output //
                  int * restrict path
)
// SYNOPSIS:                                                             
//   Log-space implementation of the Viterbi algorithm.                  
//                                                                       
// NUMERIC STABILITY:                                                    
//   This implementation is not NA-robust. Overflow and underflow are    
//   unlikely in log space and are not handled for that reason.          
//                                                                       
// ARGUMENTS:                                                            
//   'm': the number of states                                           
//   'n': the length of the sequence of observations                     
//   'log_Q': (m,m) log transition matrix.                               
//   'log_i': (m) log initial probabilities                              
//   'log_p': (m,n) log emission probabilities                           
//   'path': (n) Viterbi path.                                           
//                                                                       
// RETURN:                                                               
//   The Viterbi path.                                                   
//                                                                       
// SIDE EFFECTS:                                                         
//   Updates 'path' in place.                                            
{

   int i;        // State index.
   int j;        // State index.
   int k;        // Position of the time series.

   double thismax;
   double tmp;

   long double *array = malloc(2*m * sizeof(long double));
   int *argmax = malloc(m*n * sizeof(int));
   if (argmax == NULL || array == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      // Set path to -1 and return.
      memset(path, -1, n*sizeof(int));
      return;
   }
   long double *oldmax = array;
   long double *newmax = array + m;

   // Initial step of the algorithm.
   for (j = 0 ; j < m ; j++) newmax[j] = log_i[j+0*m] + log_p[j+0*m];
   for (k = 1 ; k < n ; k++) {
      // Set newmax to oldmax (by swapping).
      long double *swp = oldmax; oldmax = newmax; newmax = swp;
      // Viterbi recursion.
      for (j = 0 ; j < m ; j++) {
         thismax = oldmax[0] + log_Q[0+j*m];
         argmax[j+k*m] = 0;
         for (i = 1 ; i < m ; i++) {
            tmp = oldmax[i] + log_Q[i+j*m];
            if (tmp > thismax) {
               thismax = tmp;
               argmax[j+k*m] = i;
            }
         }
         newmax[j] = thismax + log_p[j+k*m];
      }
   }

   // Get final state.
   int final_state = 0;
   for (j = 1 ; j < m ; j++) if (newmax[j] > newmax[0]) final_state = j;
   path[n-1] = final_state;
   // Trace back the Viterbi path.
   for (k = n-2 ; k >= 0 ; k--) path[k] = argmax[path[k+1]+(k+1)*m];

   free(array);
   free(argmax);

   return;

}


double
block_fwdb(
   // input //
         unsigned int            m,
         unsigned int            nblocks,
   const unsigned int *          size,
   // params //
         double       * restrict Q,
         double       * restrict init,
   // output //
         double       * restrict prob,
         double       * restrict phi,
         double       * restrict sumtrans
)
// SYNOPSIS:                                                             
//   Wrapper for 'fwdb' which separates independent fragments of a       
//   time series.                                                        
//                                                                       
// ARGUMENTS:                                                            
//   'm': the number of states                                           
//   'nblocks': the number of fragments in the time series               
//   'size': (nblocks) the lengths of the fragments of the time series   
//   'Q': (m,m) transition matrix ('Q[i+j*m]' is a ij transtition)       
//   'init': (m) initial probabilities                                   
//   'prob': (n) the emssion probabilities                               
//   'phi': (m,n) probabilities given observations                       
//   'sumtrans': (m,m) sum of conditional transitions probabilties       
//                                                                       
// RETURN:                                                               
//   'void'                                                              
//                                                                       
// SIDE EFFECTS:                                                         
//   Updates 'prob', 'phi', 'sumtrans' and 'loglik' in place.            
{

   // Initialization.
   double loglik = 0.0;
   int offset = 0;
   memset(sumtrans, 0.0, m*m * sizeof(double));

   // Cycle over fragments of the time series.
   offset = 0;
   double *T = malloc(m*m * sizeof(double));
   if (T == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      return -1.0/0.0;
   }
   for (int i = 0 ; i < nblocks ; i++) {
      // NOTE: the call to `fwdb` replaces the values of 'prob' by
      // the normalized alphas.
      loglik += fwdb(m, size[i], Q, init, prob+offset, phi+offset, T);
      for (int j = 0 ; j < m*m ; j++) {
         sumtrans[j] += T[j];
      }
      offset += m * size[i];
   }

   free(T);

   return loglik;

}


int
is_undefined(
   const double * slice,
         int      m
)
// SYNOPSIS:                                                             
//   Helper function for `block_viterbi`. A set of 'm' emission
//   probabilities is considered undefined if one of them is NA,
//   or if they are all equal to -inf in log space.
{
   int n_inf = 0;
   for (int i = 0 ; i < m ; i++) {
      if (slice[i] != slice[i]) return 1;
      if (slice[i] == -INFINITY) n_inf++;
   }
   return (n_inf == m);
}


int
block_viterbi
(
   // input //
         unsigned int            m,
         unsigned int            nblocks,
   const unsigned int *          size,
   const double       * restrict Q,
   const double       * restrict init,
   const double       * restrict prob,
   // output //
                  int * restrict path
)
// SYNOPSIS:                                                             
//   Viterbi algorithm for fragmented time series. The arguments can be  
//   passed in linear or in log space.                                   
//                                                                       
// NUMERIC STABLITY:                                                     
//   This implementation is NA-robust by omission. If NAs are present    
//   at a given step, all the emission probabilities of that step are    
//   set to 0, so they do not contribute to the Viterbi path.            
//                                                                       
// ARGUMENTS:                                                            
//   'm': the number of states                                           
//   'nblocks': the number of fragments in the time series               
//   'size': (nblocks) the lengths of the fragments of the time series   
//   'Q': (m,m) transition matrix                                        
//   'init': (m) initial probabilities                                   
//   'prob': (n) emssion probabilities                                   
//   'path': (n) Viterbi path                                            
//                                                                       
// RETURN:                                                               
//   'void'                                                              
//                                                                       
// SIDE EFFECTS:                                                         
//   Updates 'path' in place.                                            
//                                                                       
// XXX NOTE XXX:                                                         
//   By using the index, the computation of the log probabilities can    
//   be made much faster. So far, this was not needed because the        
//   Viterbi algorithm is run only once per time series. In case of      
//   extremely long time series, this could make a significant           
//   difference.                                                         
{

   int n = 0;
   for (int i = 0 ; i < nblocks ; i++) n += size[i];

   double *log_Q = malloc(m*m * sizeof(double));
   double *log_i = malloc(m * sizeof(double));
   if (log_Q == NULL || log_i == NULL) {
      fprintf(stderr, "memory error: %s:%d\n", __FILE__, __LINE__);
      return 1;
   }

   // Check whether arguments are passed in linear or in log space.
   // Scan 'Q' until the first non 0 value.
   size_t idx;
   for (idx = 0 ; idx < m*m && Q[idx] == 0; idx++);
   int args_in_lin_space = Q[idx] > 0;

   // Check arguments.
   double sum_i = 0.0;
   for (size_t i = 0 ; i < m ; i++) {

      // Presence of NAs in 'init'.
      if (init[i] != init[i]) {
         fprintf(stderr, "invalid 'init' argument in '%s()'\n", __func__);
         free(log_Q);
         free(log_i);
         return -1;
      }
      // Log/lin consistency of 'init'.
      if ((init[i] >= 0) ^ args_in_lin_space) {
         fprintf(stderr, "mixed log/lin arguments in '%s()'\n", __func__);
         free(log_Q);
         free(log_i);
         return -1;
      }

      sum_i += args_in_lin_space ? init[i] : exp(init[i]);
      double sum_Q = 0.0;

      for (size_t j = 0 ; j < m ; j++) {

         // Presence of NAs in 'Q'.
         if (Q[i+j*m] != Q[i+j*m]) {
            fprintf(stderr, "invalid 'Q' argument in '%s()'\n", __func__);
            free(log_Q);
            free(log_i);
            return -1;
         }
         // Log/lin consistency of 'Q'.
         if ((Q[i+j*m] >= 0) ^ args_in_lin_space) {
            fprintf(stderr, "mixed log/lin arguments in '%s()'\n",
                  __func__);
            free(log_Q);
            free(log_i);
            return -1;
         }
         sum_Q += args_in_lin_space ? Q[i+j*m] : exp(Q[i+j*m]);
      }

      // Check that rows of 'Q' sum to 1.0.
      if (fabs(sum_Q - 1.0) > 1e-6) {
         fprintf(stderr, "'Q' is not stochastic in '%s()'\n", __func__);
         free(log_Q);
         free(log_i);
         return -1;
      }

   }

   // Check that 'init' sums to 1.0.
   if (fabs(sum_i - 1.0) > 1e-6) {
      fprintf(stderr, "'init' is not a probability in '%s()'\n", __func__);
      free(log_Q);
      free(log_i);
      return -1;
   }

   // Either way we make a copy of the emission probabilities
   // because we will replace undefined emissions by 0.0. Copying
   // 'init' and 'Q' is simpler for consistency.
   double *log_p = malloc(n*m *sizeof(double));
   if (log_p == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      free(log_Q);
      free(log_i);
      return 1;
   }
   if (args_in_lin_space) {
      // Copy variables.
      for (int i = 0 ; i < n*m ; i++) log_p[i] = log(prob[i]);
      for (int i = 0 ; i < m*m ; i++) log_Q[i] = log(Q[i]);
      for (int i = 0 ; i < m ; i++)   log_i[i] = log(init[i]);
   }
   else {
      for (int i = 0 ; i < n*m ; i++) log_p[i] = prob[i];
      for (int i = 0 ; i < m*m ; i++) log_Q[i] = Q[i];
      for (int i = 0 ; i < m ; i++)   log_i[i] = init[i];
   }

   // If an emssion probability is not available at some step, all
   // the log values are set to 0.
   int offset = 0;
   for (int i = 0 ; i < nblocks ; i++) {
      for (int k = 0 ; k < size[i] ; k++) {
         if (is_undefined(log_p + offset + k*m, m)) {
            memset(log_p + offset + k*m, 0, m * sizeof(double));
         }
      }
      offset += m * size[i];
   }

   // NOTE: the offset is not the same in 'path' and 'log_p' because
   // of their dimensions (explains 'm*offset' in the case of 'log_p').
   offset = 0;
   for (int i = 0 ; i < nblocks ; i++) {
      viterbi(m, size[i], log_Q, log_i, log_p+m*offset, path+offset);
      offset += size[i];
   }

   free(log_Q);
   free(log_i);
   free(log_p);

   return 0;

}
