#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "../../src/jahmm.h"


// Register the method 'jahmm_R_call()'.
SEXP jahmm_R_call (SEXP);

R_CallMethodDef callMethods[] = {
   {"jahmm_R_call", (DL_FUNC) &jahmm_R_call, 1},
   {NULL, NULL, 0},
};

void R_init_jahmm(DllInfo *info) {
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

SEXP
jahmm_R_call
(
   SEXP Y
)
{

   const unsigned int r = length(Y) - 1;
   const unsigned int n = length(VECTOR_ELT(Y, 0));

   // Get the values.
   int *y = malloc(n*r * sizeof(int));
   if (y == NULL) {
      Rprintf("Rjahmm memory error %s:%d\n", __FILE__, __LINE__);
      return R_NilValue;
   }

   for (size_t i = 1 ; i < r+1 ; i++) {
      int *v = INTEGER(coerceVector(VECTOR_ELT(Y, i), INTSXP));
      for (size_t j = 0 ; j < n ; j++) y[i-1+j*r] = v[j];
   }

   // Get the blocks.
   histo_t *histo = new_histo();
   if (histo == NULL) {
      Rprintf("Rjahmm memory error %s:%d\n", __FILE__, __LINE__);
      return R_NilValue;
   }

   int id_array[1] = {-1};
   int *id = id_array;
   char prev[256] = {0};

   SEXP S = VECTOR_ELT(Y, 0);
   switch (TYPEOF(S)) {
      case STRSXP:
         for (size_t i = 0 ; i < n ; i++) {
            const char *block = CHAR(STRING_ELT(S, i));
            if (strcmp(prev, block) != 0) {
               strncpy(prev, block, 256);
               (*id)++;
            }
            histo_push(&histo, *id);
         }
         break;
      case INTSXP:
         id = INTEGER(S);
         for (size_t i = 0 ; i < n ; i++) {
            histo_push(&histo, id[i]);
         }
   }

   tab_t *tab = compress_histo(histo);
   if (tab == NULL) {
      Rprintf("Rjahmm memory error %s:%d\n", __FILE__, __LINE__);
      return R_NilValue;
   }
   ChIP_t *ChIP = new_ChIP(r, tab->size, y, tab->num);
   free(histo);
   free(tab);

   jahmm_t * jahmm = do_jahmm(ChIP);

   if (jahmm == NULL) {
      Rprintf("Rjahmm error\n");
      return R_NilValue;
   }

   // Jahmm uses 3 states.
   const unsigned int m = 3;

   SEXP Q;
   PROTECT(Q = allocVector(REALSXP, m*m));
   memcpy(REAL(Q), jahmm->Q, m*m * sizeof(double));

   SEXP A;
   PROTECT(A = allocVector(REALSXP, 1));
   *REAL(A) =  jahmm->a;

   SEXP PI_;
   PROTECT(PI_ = allocVector(REALSXP, 1));
   *REAL(PI_) =  jahmm->pi;

   // The parameters 'p' are coded "row-wise". Since R
   // objects are coded "column-wise" we neeed to
   // disentangle them.
   SEXP P;
   PROTECT(P = allocVector(REALSXP, m*(r+1)));
   for (size_t i = 0 ; i < r+1 ; i++) {
   for (size_t j = 0 ; j < m ; j++) {
      REAL(P)[j+i*m] = jahmm->p[i+j*(r+1)];
   }
   }

   // For reason of cache-friendlyness, 'phi' and 'pem'
   // are also coded "row-wise".
   SEXP PHI;
   SEXP PEM;
   PROTECT(PHI = allocVector(REALSXP, m*n));
   PROTECT(PEM = allocVector(REALSXP, m*n));
   for (size_t i = 0 ; i < n ; i++) {
   for (size_t j = 0 ; j < m ; j++) {
      REAL(PEM)[i+j*n] = jahmm->pem[j+i*m];
      REAL(PHI)[i+j*n] = jahmm->phi[j+i*m];
   }
   }

   SEXP PATH;
   PROTECT(PATH = allocVector(INTSXP, n));
   memcpy(INTEGER(PATH), jahmm->path, n * sizeof(int));

   SEXP L;
   PROTECT(L = allocVector(REALSXP, 1));
   *REAL(L) = jahmm->l;

   SEXP dimQ;
   PROTECT(dimQ = allocVector(INTSXP, 2));
   INTEGER(dimQ)[0] = m;
   INTEGER(dimQ)[1] = m;
   setAttrib(Q, R_DimSymbol, dimQ);

   SEXP dimP;
   PROTECT(dimP = allocVector(INTSXP, 2));
   INTEGER(dimP)[0] = m;
   INTEGER(dimP)[1] = r+1;
   setAttrib(P, R_DimSymbol, dimP);

   SEXP dimPHI;
   PROTECT(dimPHI = allocVector(INTSXP, 2));
   INTEGER(dimPHI)[0] = n;
   INTEGER(dimPHI)[1] = m;
   setAttrib(PHI, R_DimSymbol, dimPHI);

   SEXP dimPEM;
   PROTECT(dimPEM = allocVector(INTSXP, 2));
   INTEGER(dimPEM)[0] = n;
   INTEGER(dimPEM)[1] = m;
   setAttrib(PEM, R_DimSymbol, dimPEM);

   SEXP RETLIST;
   PROTECT(RETLIST = allocVector(VECSXP, 8));
   SET_VECTOR_ELT(RETLIST, 0, Q);
   SET_VECTOR_ELT(RETLIST, 1, A);
   SET_VECTOR_ELT(RETLIST, 2, PI_);
   SET_VECTOR_ELT(RETLIST, 3, P);
   SET_VECTOR_ELT(RETLIST, 4, PHI);
   SET_VECTOR_ELT(RETLIST, 5, PEM);
   SET_VECTOR_ELT(RETLIST, 6, PATH);
   SET_VECTOR_ELT(RETLIST, 7, L);

   UNPROTECT(13);

   destroy_jahmm_all(jahmm);

   return RETLIST;

}
