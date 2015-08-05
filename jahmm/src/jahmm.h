#include <stdint.h>
#include "hmm.h"
#include "utils.h"
#include "zinb.h"

#ifndef _JAHMM_HEADER
#define _JAHMM_HEADER

#define BW_MAXITER 100     // BW iterations //
#define BT_MAXITER 20      // Backtrack iterations //
#define TOLERANCE 1e-6

struct ChIP_t;
struct jahmm_t;

typedef struct ChIP_t ChIP_t;
typedef struct jahmm_t jahmm_t;


struct ChIP_t {
   size_t          r;      // number of dimensions of 'y' //
   unsigned int    nb;     // number of blocks //
   int           * y;      // observations //
   unsigned int    size[]; // size of the blocks //
};

struct jahmm_t {
   unsigned int    m;      // number of states //
   ChIP_t        * ChIP;   // observations //
   double        * Q;      // transitions //
   double          a;      // emission par //
   double          pi;     // emission par //
   double        * p;      // emission par //
   double        * phi;    // posterior probs //
   double        * pem;    // emission probs //
   double          l;      // log-likelihood //
   int           * path;   // Viterbi path //
   int             iter;   // number of BW iterations //
};

void      bw_zinm(jahmm_t *);
void      destroy_jahmm_all(jahmm_t *);
jahmm_t * do_jahmm(ChIP_t *);
ChIP_t  * new_ChIP(unsigned int, unsigned int, int *,
            const unsigned int *);
jahmm_t * new_jahmm(unsigned int, ChIP_t *);
unsigned int nobs(const ChIP_t *);
ChIP_t  * read_file(FILE *);
void      set_jahmm_par(jahmm_t *, const double *, double,
            double, const double *);
void      update_trans(size_t, double *, const double *);
void      zinm_prob(jahmm_t *, const int *, int, double *);
#endif
