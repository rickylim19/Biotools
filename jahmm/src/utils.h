#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "xxhash.h"

#ifndef _ZINM_UTILS_HEADER
#define _ZINM_UTILS_HEADER

#define HISTO_INIT_SIZE 128

struct histo_t;
struct tab_t;

typedef struct histo_t histo_t;
typedef struct tab_t tab_t;

struct histo_t {
   size_t size;
   unsigned int num[];
};

struct tab_t {
   size_t size;
   unsigned int *num;
   unsigned int val[];
};

int       indexts (int, int, const int *, int *);
tab_t   * compress_histo(histo_t *);
double    digamma(double);
int       histo_push(histo_t **, size_t);
tab_t   * init_tab(size_t, unsigned int *, unsigned int *);
histo_t * new_histo(void);
void      rowsums(int *, unsigned int, unsigned int, int*);
double    trigamma(double);
tab_t   * tabulate(int *, unsigned int);

#endif
