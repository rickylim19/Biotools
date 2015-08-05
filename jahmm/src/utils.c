#include "utils.h"

#define U32 uint32_t

// Globals variable.
U32 *xxhash;

void
rowsums
(
   int *x,
   unsigned int dim,
   unsigned int nobs,
   int *sums
)
{

   // Assume data is given by rows in x.
   // Negative values are considered NA.
   memset(sums, 0, nobs * sizeof(int));
   for (size_t i = 0 ; i < nobs ; i++) {
   for (size_t j = 0 ; j < dim ; j++) {
      if (x[i+j*dim] < 0) {
         sums[i] = -1;
         break;
      }
      sums[i] += x[j+i*dim];
   }
   }
}


tab_t *
tabulate
(
   int *x,
   unsigned int nobs
)
{

   histo_t *histo = new_histo();
   if (histo == NULL) return NULL;
   for (size_t i = 0 ; i < nobs ; i++) {
      // Skip negative values (used for NA).
      if (x[i] < 0) continue;
      if (histo_push(&histo, x[i])) {
         free(histo);
         return NULL;
      }
   }

   tab_t *tab = compress_histo(histo);
   free(histo);
   return tab;

}


histo_t *
new_histo
(void)
{

   size_t initsize = sizeof(histo_t) + HISTO_INIT_SIZE * sizeof(int);
   histo_t *new = calloc(1, initsize);
   if (new == NULL) {
      fprintf(stderr, "memory error: %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }
   new->size = HISTO_INIT_SIZE;

   return new;

}


int
histo_push
(
   histo_t **histo_addr,
   size_t val
)
{

   // Convenience variable.
   histo_t *histo = *histo_addr;
   if (val >= histo->size) {
      size_t newsize = 2*val * (sizeof(int));
      histo_t *new = realloc(histo, sizeof(histo_t) + newsize);
      if (new == NULL) {
         fprintf(stderr, "memory error: %s:%d\n", __FILE__, __LINE__);
         return 1;
      }
      *histo_addr = histo = new;
      size_t added_size = (2*val - histo->size) * sizeof(int);
      memset(histo->num + histo->size, 0, added_size);
      histo->size = 2*val;
   }

   histo->num[val]++;
   return 0;

}


tab_t *
init_tab
(
   size_t size,
   unsigned int *val,
   unsigned int *num
)
{
   size_t extra = 2*size * sizeof(int);
   tab_t *new = malloc(sizeof(tab_t) + extra);
   if (new == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }

   new->size = size;
   new->num = new->val + size;
   memcpy(new->val, val, size * sizeof(int));
   memcpy(new->num, num, size * sizeof(int));

   return new;

}


tab_t *
compress_histo
(
   histo_t *histo
)
{

   size_t size = 0;
   for (size_t i = 0 ; i < histo->size ; i++) {
      size += (histo->num[i] > 0);
   }
   tab_t *new = malloc(sizeof(tab_t) + 2*size * sizeof(int));
   if (new == NULL) {
      fprintf(stderr, "memory error: %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }
   new->size = size;
   new->num = new->val + size;

   size_t j = 0;
   for (size_t i = 0 ; i < histo->size ; i++) {
      if (histo->num[i] > 0) {
         new->val[j] = i;
         new->num[j] = histo->num[i];
         j++;
      }
   }

   return new;

}


int
stblcmp
(
   const void *a,
   const void *b
)
// SYNOPSIS:                                                             
//   Comparison function for stable sort on hashes. Used to sort         
//   addresses of global pointer 'xxhash'.                               
{
   U32 A = xxhash[*(int *)a];
   U32 B = xxhash[*(int *)b];
   if (A > B) return 1;
   if (A < B) return -1;
   // Hashes are identical, compare addresses for stable sort.
   if (*(int *)a > *(int *)b)   return  1;
   return -1;
}


int
indexts
(
         int n,
         int r,
   const int *ts,
         int *index
)
// SYNOPSIS:
// Index the time series using xxhash and return the index of the
// first all-0 observation.
{
   int i;
   xxhash = malloc(n*sizeof(U32));
   int *addr = malloc(n * sizeof(int));
   if (xxhash == NULL || addr == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      return -1;
   }

   for (i = 0 ; i < n ; i++) addr[i] = i;

   // Compute xxhash digests.
   for (i = 0 ; i < n ; i++) {
      xxhash[i] = XXH32(ts+i*r, r*sizeof(int), 0);
   }

   // Compute the xxhash digest of all 0s. We will need it
   // to return the index of the small such observation.
   int *all0 = calloc(r, sizeof(int));
   if (all0 == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      return -1;
   }
   U32 xxhash0 = XXH32(all0, r*sizeof(int), 0);
   free(all0);

   // Stable sort array indices on digests order. Stability is
   // important because we want the first occurrence of an
   // observation to point to its own index, and all subsequent
   // occurrences to point to it as well. If the reference index
   // was not the first occurence in array order, it would cause
   // difficulties for the purpose of computing emission
   // probabilities.
   qsort(addr, n, sizeof(int), stblcmp);

   int current = 0;
   int index0 = -1;
   index[addr[0]] = addr[0];
   if (xxhash[addr[0]] == xxhash0) index0 = addr[0];
   for (i = 1 ; i < n ; i++) {
      if (xxhash[addr[i]] == xxhash[current]) {
         index[addr[i]] = current;
      }
      else {
         current = index[addr[i]] = addr[i];
         if (xxhash[addr[i]] == xxhash0) index0 = addr[i];
      }
   }

   free(xxhash);
   free(addr);

   return index0;

}


// The code below was copied from the following link
// http://pmtksupport.googlecode.com/svn/trunk/lightspeed2.3/util.c
// Written by Tom Minka (unless otherwise noted).

/* The digamma function is the derivative of gammaln.

   Reference:
    J Bernardo,
    Psi ( Digamma ) Function,
    Algorithm AS 103,
    Applied Statistics,
    Volume 25, Number 3, pages 315-317, 1976.

    From http://www.psc.edu/~burkardt/src/dirichlet/dirichlet.f
    (with modifications for negative numbers and extra precision)
*/


double
digamma
(
   double x
)
{
  double neginf = -INFINITY;
  static const double c = 12,
    digamma1 = -0.57721566490153286,
    trigamma1 = 1.6449340668482264365, /* pi^2/6 */
    s = 1e-6,
    s3 = 1./12,
    s4 = 1./120,
    s5 = 1./252,
    s6 = 1./240,
    s7 = 1./132;
    //s8 = 691./32760,
    //s9 = 1./12,
    //s10 = 3617./8160;
  double result;
  /* Illegal arguments */
  if((x == neginf) || isnan(x)) {
    return NAN;
  }
  /* Singularities */
  if((x <= 0) && (floor(x) == x)) {
    return neginf;
  }
  /* Negative values */
  /* Use the reflection formula (Jeffrey 11.1.6):
   * digamma(-x) = digamma(x+1) + pi*cot(pi*x)
   *
   * This is related to the identity
   * digamma(-x) = digamma(x+1) - digamma(z) + digamma(1-z)
   * where z is the fractional part of x
   * For example:
   * digamma(-3.1) = 1/3.1 + 1/2.1 + 1/1.1 + 1/0.1 + digamma(1-0.1)
   *               = digamma(4.1) - digamma(0.1) + digamma(1-0.1)
   * Then we use
   * digamma(1-z) - digamma(z) = pi*cot(pi*z)
   */
  if(x < 0) {
    return digamma(1-x) + M_PI / tan(-M_PI*x);
  }
  /* Use Taylor series if argument <= S */
  if(x <= s) return digamma1 - 1/x + trigamma1*x;
  /* Reduce to digamma(X + N) where (X + N) >= C */
  result = 0;
  while(x < c) {
    result -= 1/x;
    x++;
  }
  /* Use de Moivre's expansion if argument >= C */
  /* This expansion can be computed in Maple via asympt(Psi(x),x) */
  if(x >= c) {
    double r = 1/x;
    result += log(x) - 0.5*r;
    r *= r;
    result -= r * (s3 - r * (s4 - r * (s5 - r * (s6 - r * s7))));
  }
  return result;
}

/* The trigamma function is the derivative of the digamma function.

   Reference:

    B Schneider,
    Trigamma Function,
    Algorithm AS 121,
    Applied Statistics, 
    Volume 27, Number 1, page 97-99, 1978.

    From http://www.psc.edu/~burkardt/src/dirichlet/dirichlet.f
    (with modification for negative arguments and extra precision)
*/


double trigamma(double x)
{
  double neginf = -INFINITY,
    small = 1e-4,
    large = 8,
    trigamma1 = 1.6449340668482264365, /* pi^2/6 = Zeta(2) */
    tetragamma1 = -2.404113806319188570799476,  /* -2 Zeta(3) */
    b2 =  1./6,  /* B_2 */
    b4 = -1./30, /* B_4 */
    b6 =  1./42, /* B_6 */
    b8 = -1./30, /* B_8 */
    b10 = 5./66; /* B_10 */
  double result;
  /* Illegal arguments */
  if((x == neginf) || isnan(x)) {
    return NAN;
  }
  /* Singularities */
  if((x <= 0) && (floor(x) == x)) {
    return neginf;
  }
  /* Negative values */
  /* Use the derivative of the digamma reflection formula:
   * -trigamma(-x) = trigamma(x+1) - (pi*csc(pi*x))^2
   */
  if(x < 0) {
    result = M_PI/sin(-M_PI*x);
    return -trigamma(1-x) + result*result;
  }
  /* Use Taylor series if argument <= small */
  if(x <= small) {
    return 1/(x*x) + trigamma1 + tetragamma1*x;
  }
  result = 0;
  /* Reduce to trigamma(x+n) where ( X + N ) >= B */
  while(x < large) {
    result += 1/(x*x);
    x++;
  }
  /* Apply asymptotic formula when X >= B */
  /* This expansion can be computed in Maple via asympt(Psi(1,x),x) */
  if(x >= large) {
    double r = 1/(x*x);
    result += 0.5*r + (1 + r*(b2 + r*(b4 + r*(b6 + r*(b8 + r*b10)))))/x;
  }
  return result;
}

/*
double
mean
(
   const int *yz,
   int n,
   int r
)
// SYNOPSIS:                                                             
//   Compute the mean of the first column of an integer array. To com-   
//   pute the mean of another column, call as `mean(yz+1, n, r)`, where  
//   1 is for column 2 etc.                                              
//                                                                       
// PARAMETERS:                                                           
//   'yz': (r,n) the integer array                                       
//   'n': row number of the array                                        
//   'r': column number of the array                                     
//                                                                       
// RETURN:                                                               
//   The mean as a 'double'.                                             
{
   double sum = 0.0;
   int n_obs_no_NA = 0;

   for (int k = 0 ; k < n ; k++) {
      // Casting NA to integer gives -2147483648, which is the 
      // largest negative value stored in 'int'. Here I test for
      // NA by wrapping around.
      if (yz[r*k] == INT_MIN) continue;
      sum += yz[r*k];
      n_obs_no_NA++;
   }
   // The result can be 0/0, which is 'nan'.
   return sum / n_obs_no_NA;
}
*/

/*
int *
histsum
(
   const int *yz,
   int n,
   int r
)
// SYNOPSIS:                                                             
//   Compute the histogram of the row-wise sum of an integer array.      
//                                                                       
// INPUT:                                                                
//   The presence of a negative value makes the whole row ignored.       
//                                                                       
// PARAMETERS:                                                           
//   'yz': (r,n) the integer array                                       
//   'n': row number of the array                                        
//   'r': column number of the array                                     
//                                                                       
// RETURN:                                                               
//   A pointer of 'int' to the histogram.                                
{
   int size = 1024;
   int *counts = calloc(size, sizeof(int));
   if (counts == NULL) {
      fprintf(stderr, "memory error (hist 1)\n");
      return NULL;
   }

   int maxval = 0;
   for (int k = 0 ; k < n ; k++) {
      int sum = 0;
      for (int i = 0 ; i < r ; i++){
         if (yz[i+k*r] < 0) {
            sum = -1;
            break;
         }
         sum += yz[i+k*r];
      }
      if (sum < 0) continue;
      if (sum > maxval) maxval = sum;
      if (sum > size-2) {
         int newsize;
         for (newsize = size ; sum > newsize-2 ; newsize *= 2);
         int *newcounts = realloc(counts, newsize * sizeof(int));
         if (newcounts == NULL) {
            fprintf(stderr, "memory error (hist 2)\n");
            return NULL;
         }
         else {
            // Extra memory must be initialized to 0.
            counts = newcounts;
            for (int j = size ; j < newsize ; j++) counts[j] = 0;
            size = newsize;
         }
      }
      counts[sum]++;
   }
   // Add the sentinel.
   counts[maxval+1] = -1;
   return counts;

}
*/


