#include "hmm.h"
#include "jahmm.h"
#include "unittest.h"
#include "utils.h"
#include "zinb.h"

void
test_new_histo
(void)
{

   histo_t *histo;
   histo = new_histo();
   test_assert_critical(histo != NULL);
   for (size_t i = 0 ; i < HISTO_INIT_SIZE ; i++) {
      test_assert(histo->num[i] == 0);
   }
   free(histo);

   redirect_stderr();
   set_alloc_failure_rate_to(1.0);
   histo = new_histo();
   reset_alloc();
   unredirect_stderr();
   test_assert(histo == NULL);
   test_assert(strcmp(caught_in_stderr(), "") != 0);

   return;

}

void
test_histo_push
(void)
{

   histo_t *histo;
   histo = new_histo();
   test_assert_critical(histo != NULL);
   for (size_t i = 0 ; i < (1 << 16) ; i++) {
      test_assert(histo_push(&histo, i) == 0);
   }
   test_assert(histo->size == (1 << 16));
   for (size_t i = 0 ; i < (1 << 16) ; i++) {
      test_assert(histo->num[i] == 1);
   }
   free(histo);

   histo = new_histo();
   test_assert_critical(histo != NULL);
   for (size_t i = (1 << 16) ; i > 0 ; i--) {
      test_assert(histo_push(&histo, i-1) == 0);
   }
   test_assert(histo->size == 2*((1 << 16)-1));
   for (size_t i = 0 ; i < (1 << 16) ; i++) {
      test_assert(histo->num[i] == 1);
   }
   free(histo);

   histo = new_histo();
   test_assert_critical(histo != NULL);
   redirect_stderr();
   set_alloc_failure_rate_to(1.0);
   test_assert(histo_push(&histo, HISTO_INIT_SIZE) == 1);
   reset_alloc();
   unredirect_stderr();
   test_assert(strcmp(caught_in_stderr(), "") != 0);
   free(histo);

   return;

}

void
test_compress_histo
(void)
{

   histo_t *histo = new_histo();
   test_assert_critical(histo != NULL);
   for (size_t i = 0 ; i < 4096 ; i += 2) {
      test_assert(histo_push(&histo, i) == 0);
   }

   tab_t *tab;
   tab = compress_histo(histo);
   test_assert_critical(tab != NULL);
   test_assert(tab->size == 2048);
   for (size_t i = 0 ; i < tab->size ; i++) {
      test_assert(tab->val[i] == 2*i);
      test_assert(tab->num[i] == 1);
   }

   free(tab);

   redirect_stderr();
   set_alloc_failure_rate_to(1.0);
   tab = compress_histo(histo);
   reset_alloc();
   unredirect_stderr();
   test_assert(tab == NULL);
   test_assert(strcmp(caught_in_stderr(), "") != 0);

   free(histo);
   return;

}

void
test_tabulate
(void)
{

   int x1[12] = {1,2,3,4,1,2,3,4,1,2,3,4};
   tab_t *tab;
   tab = tabulate(x1, 12);
   test_assert_critical(tab != NULL);
   test_assert(tab->size == 4);
   test_assert(tab->val[0] == 1);
   test_assert(tab->val[1] == 2);
   test_assert(tab->val[2] == 3);
   test_assert(tab->val[3] == 4);
   for (int i = 0 ; i < 4 ; i++) {
      test_assert(tab->num[i] == 3);
   }
   free(tab);

   int x2[4] = {4096,2048,1024,512};
   tab = tabulate(x2, 4);
   test_assert_critical(tab != NULL);
   test_assert(tab->size == 4);
   test_assert(tab->val[0] == 512);
   test_assert(tab->val[1] == 1024);
   test_assert(tab->val[2] == 2048);
   test_assert(tab->val[3] == 4096);
   for (int i = 0 ; i < 4 ; i++) {
      test_assert(tab->num[i] == 1);
   }
   free(tab);

   int x3[5] = {-1,2,3,-1,2};
   tab = tabulate(x3, 5);
   test_assert_critical(tab != NULL);
   test_assert(tab->size == 2);
   test_assert(tab->val[0] == 2);
   test_assert(tab->val[1] == 3);
   test_assert(tab->num[0] == 2);
   test_assert(tab->num[1] == 1);
   free(tab);

   return;

}


void
test_eval_nb_f
(void)
{

   // 0:14, 1:5, 2:4, 3:1, 5:1
   int x[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                     1,1,1,1,1,2,2,2,2,3,5 };
   tab_t *tab = tabulate(x, 25);
   test_assert(fabs(eval_nb_f(1.0, tab)+0.12747262) < 1e-6);
   test_assert(fabs(eval_nb_f(1.1, tab)+0.24215981) < 1e-6);
   test_assert(fabs(eval_nb_f(1.2, tab)+0.31636395) < 1e-6);
   test_assert(fabs(eval_nb_f(1.3, tab)+0.36350700) < 1e-6);
   test_assert(fabs(eval_nb_f(1.4, tab)+0.39225466) < 1e-6);
   test_assert(fabs(eval_nb_f(1.5, tab)+0.40834322) < 1e-6);
   test_assert(fabs(eval_nb_f(2.0, tab)+0.39975512) < 1e-6);

   free(tab);

   return;

}


void
test_eval_nb_dfda
(void)
{

   // 0:14, 1:5, 2:4, 3:1, 5:1
   int x[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                    1,1,1,1,1,2,2,2,2,3,5 };
   tab_t *tab = tabulate(x, 25);
   test_assert(fabs(eval_nb_dfda(1.0, tab)+1.41167874) < 1e-6);
   test_assert(fabs(eval_nb_dfda(1.2, tab)+0.58911102) < 1e-6);
   test_assert(fabs(eval_nb_dfda(1.3, tab)+0.36790287) < 1e-6);
   test_assert(fabs(eval_nb_dfda(1.4, tab)+0.21643981) < 1e-6);
   test_assert(fabs(eval_nb_dfda(1.5, tab)+0.11168877) < 1e-6);
   test_assert(fabs(eval_nb_dfda(2.0, tab)-0.08773865) < 1e-6);

   free(tab);

   return;

}


void
test_eval_zinb_f
(void)
{

   test_assert(fabs(eval_zinb_f(1.0, .5, 1, 2.0)) < 1e-6);
   test_assert(fabs(eval_zinb_f(1.0, .5, 2, 4.0)) < 1e-6);
   test_assert(fabs(eval_zinb_f(1.3, .7, 141, 8.3)-678.0838351) < 1e-6);

   return;

}


void
test_eval_zinb_g
(void)
{

   tab_t *tab;
   // 0:14, 1:5, 2:4, 3:1, 5:1
   int x1[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  1,1,1,1,1,2,2,2,2,3,5 };
   tab = tabulate(x1, 25);
   test_assert_critical(tab != NULL);
   test_assert(fabs(eval_zinb_g(1, .5, tab)+0.1325713) < 1e-6);
   free(tab);

   return;

}


void
test_eval_zinb_dfda
(void)
{

   test_assert(fabs(eval_zinb_dfda(1, .5, 1)-1.2274112) < 1e-6);
   test_assert(fabs(eval_zinb_dfda(1, .5, 9)-11.0467015) < 1e-6);
   test_assert(fabs(eval_zinb_dfda(2, .5, 9)-12.9096451) < 1e-6);
   test_assert(fabs(eval_zinb_dfda(2, .3, 9)-25.1159846) < 1e-6);
   test_assert(fabs(eval_zinb_dfda(2.4, .3, 9)-26.3620864) < 1e-6);

   return;

}


void
test_eval_zinb_dfdp
(void)
{

   test_assert(eval_zinb_dfdp(1, .5, 0, 0) == 0.0);
   test_assert(eval_zinb_dfdp(1, .5, 1, 0) == 0.0);
   test_assert(fabs(eval_zinb_dfdp(2, .5, 1, 0)+3.5555555) < 1e-6);
   test_assert(fabs(eval_zinb_dfdp(2, .3, 1, 0)+19.5896899) < 1e-6);
   test_assert(fabs(eval_zinb_dfdp(2, .3, 9, 0)+176.3072092) < 1e-6);
   test_assert(fabs(eval_zinb_dfdp(2, .3, 9, 1.7)+179.7765970) < 1e-6);

   return;

}


void
test_eval_zinb_dgda
(void)
{

   tab_t *tab;
   // 0:14, 1:5, 2:4, 3:1, 5:1
   int x1[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  1,1,1,1,1,2,2,2,2,3,5 };
   tab = tabulate(x1, 25);
   test_assert_critical(tab != NULL);
   test_assert(fabs(eval_zinb_dgda(1, .5, tab)+2.2547559) < 1e-6);
   test_assert(fabs(eval_zinb_dgda(2, .5, tab)+1.2605630) < 1e-6);
   test_assert(fabs(eval_zinb_dgda(2, .3, tab)+1.8764955) < 1e-6);
   free(tab);

   return;

}


void
test_ll_zinb
(void)
{

   tab_t *tab;
   // 0:14, 1:5, 2:4, 3:1, 5:1
   int x1[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  1,1,1,1,1,2,2,2,2,3,5 };
   tab = tabulate(x1, 25);
   test_assert_critical(tab != NULL);
   test_assert(fabs(ll_zinb(1, .5, 1, tab)+22.5329303) < 1e-6);
   test_assert(fabs(ll_zinb(1, .5, .7, tab)+22.7832550) < 1e-6);
   test_assert(fabs(ll_zinb(2, .5, .7, tab)+23.7608409) < 1e-6);
   test_assert(fabs(ll_zinb(2, .3, .7, tab)+31.6978553) < 1e-6);
   free(tab);

   return;

}


void
test_nb_est_alpha
(void)
{

   tab_t *tab;
   // 0:14, 1:5, 2:4, 3:1, 5:1
   int x1[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  1,1,1,1,1,2,2,2,2,3,5 };
   tab = tabulate(x1, 25);
   test_assert_critical(tab != NULL);
   test_assert(fabs(nb_est_alpha(tab)-0.9237) < 1e-3);
   free(tab);

   // 0:27, 1:12, 2:8, 3:1, 4:1, 5:1
   int x2[50] = {3,0,1,2,0,0,1,0,0,0,0,1,1,0,0,1,2,2,0,0,0,1,2,
      0, 0,0,0,0,4,0,0,0,1,5,1,0,1,2,1,2,2,2,0,0,0,1,0,1,0,0};
   tab = tabulate(x2, 50);
   test_assert_critical(tab != NULL);
   test_assert(fabs(nb_est_alpha(tab)-1.3436) < 1e-3);
   free(tab);

   // 0:12, 1:7, 2:13, 3:4, 4:6, 5:2, 6:1, 7:3, 8:1, 9:1
   int x3[50] = {4,5,2,1,2,4,2,2,0,4,2,1,3,6,0,0,7,3,0,8,4,2,0,
      0,2,3,2,3,7,9,2,4,0,4,2,0,0,2,5,1,1,2,1,0,0,0,1,2,1,7};
   tab = tabulate(x3, 50);
   test_assert_critical(tab != NULL);
   test_assert(fabs(nb_est_alpha(tab)-1.7969) < 1e-3);
   free(tab);

   // 0:39, 1:8, 2:2, 3:1
   int x4[50] = {1,0,0,0,0,0,3,1,0,1,0,0,0,0,0,0,0,1,0,0,0,2,0,
      2,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0};
   tab = tabulate(x4, 50);
   test_assert_critical(tab != NULL);
   test_assert(fabs(nb_est_alpha(tab)-0.7073) < 1e-3);
   free(tab);

   // 0:59, 1:83, 2:99, 3:67, 4:67, 5:49, 6:27, 7:22, 8:11, 9:6
   // 10:6, 11:3, 12:2, 13:3
   int x5[500] = {1,0,0,1,1,2,1,0,5,7,1,3,3,1,6,0,2,5,7,0,5,2,1,
       10,5,3,4,5,7,0,8,6,3,0,2,1,1,0,2,3,7,2,3,2,2,1,0,4,4,2,4,2,
       0,6,3,2,5,2,1,4,3,4,2,2,5,3,2,0,2,8,1,3,1,7,5,1,4,1,1,0,2,
       2,4,1,1,1,4,1,3,4,4,10,5,2,0,7,1,6,1,3,6,4,0,2,4,1,12,2,5,
       6,5,4,1,11,0,1,3,2,4,2,0,2,3,4,0,2,9,9,7,4,2,1,3,3,3,4,2,9,
       2,4,3,2,2,4,2,5,3,0,1,3,2,0,3,3,4,1,3,3,5,7,3,3,2,1,5,5,4,
       6,1,1,1,2,9,5,1,2,4,0,2,1,0,3,2,4,3,1,4,2,1,4,1,6,0,6,5,3,
       5,2,0,1,2,1,0,5,3,2,7,6,4,3,2,5,7,5,5,1,1,3,10,2,0,5,0,1,2,
       0,5,1,2,3,6,4,0,3,1,2,2,4,3,0,3,2,5,4,10,1,2,4,4,2,13,4,3,
       1,5,4,8,5,6,2,3,4,3,1,5,5,1,8,2,0,5,7,3,2,2,4,2,3,1,5,3,7,
       13,1,4,7,5,5,0,3,0,4,2,3,1,2,4,2,8,1,2,5,6,1,1,0,7,2,2,3,5,
       12,2,2,2,0,3,3,4,0,2,5,1,10,0,7,6,5,0,11,2,3,7,3,5,4,2,1,2,
       4,0,2,2,2,0,6,2,3,4,2,3,7,3,5,2,5,0,4,4,6,3,1,2,7,3,0,2,5,
       7,2,2,0,0,0,6,3,0,1,1,5,5,2,6,2,4,6,0,1,2,3,2,2,2,3,4,1,1,
       4,0,2,0,1,3,4,1,2,2,3,1,4,4,3,4,4,1,5,2,13,4,10,5,6,1,0,5,
       0,0,5,6,0,1,8,5,1,3,1,8,1,8,1,6,7,2,8,2,2,3,3,0,4,2,1,9,6,
       0,6,7,1,8,2,2,1,11,3,0,4,2,5,1,6,8,3,4,7,0,4,2,4,1,1,1,6,0,
       4,4,6,2,1,3,1,0,4,9,3,1,4,2,2,0,1};
   tab = tabulate(x5, 500);
   test_assert_critical(tab != NULL);
   test_assert(fabs(nb_est_alpha(tab)-3.0057) < 1e-3);
   free(tab);
}


void
test_mle_nb
(void)
{

   // These test cases have been verified with R.
   zinb_par_t *par;
   
   // 0:14, 1:5, 2:4, 3:1, 5:1
   int x1[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  1,1,1,1,1,2,2,2,2,3,5 };
   par = mle_nb(x1, 25);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->a - 0.9237) < 1e-3);
   test_assert(par->pi == 0.0);
   test_assert(fabs(par->p - 0.5237) < 1e-3);
   free(par);

   // 0:27, 1:12, 2:8, 3:1, 4:1, 5:1
   int x2[50] = {3,0,1,2,0,0,1,0,0,0,0,1,1,0,0,1,2,2,0,0,0,1,2,
      0, 0,0,0,0,4,0,0,0,1,5,1,0,1,2,1,2,2,2,0,0,0,1,0,1,0,0};
   par = mle_nb(x2, 50);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->a - 1.3436) < 1e-3);
   test_assert(par->pi == 0.0);
   test_assert(fabs(par->p - 0.6267) < 1e-3);
   free(par);

   // 0:12, 1:7, 2:13, 3:4, 4:6, 5:2, 6:1, 7:3, 8:1, 9:1
   int x3[50] = {4,5,2,1,2,4,2,2,0,4,2,1,3,6,0,0,7,3,0,8,4,2,0,
      0,2,3,2,3,7,9,2,4,0,4,2,0,0,2,5,1,1,2,1,0,0,0,1,2,1,7};
   par = mle_nb(x3, 50);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->a - 1.7969) < 1e-3);
   test_assert(par->pi == 0.0);
   test_assert(fabs(par->p - 0.4221) < 1e-3);
   free(par);

   // 0:39, 1:8, 2:2, 3:1
   int x4[50] = {1,0,0,0,0,0,3,1,0,1,0,0,0,0,0,0,0,1,0,0,0,2,0,
      2,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0};
   par = mle_nb(x4, 50);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->a - 0.7073) < 1e-3);
   test_assert(par->pi == 0.0);
   test_assert(fabs(par->p - 0.7021) < 1e-3);
   free(par);

   // 0:59, 1:83, 2:99, 3:67, 4:67, 5:49, 6:27, 7:22, 8:11, 9:6
   // 10:6, 11:3, 12:2, 13:3
   int x5[500] = {1,0,0,1,1,2,1,0,5,7,1,3,3,1,6,0,2,5,7,0,5,2,1,
       10,5,3,4,5,7,0,8,6,3,0,2,1,1,0,2,3,7,2,3,2,2,1,0,4,4,2,4,2,
       0,6,3,2,5,2,1,4,3,4,2,2,5,3,2,0,2,8,1,3,1,7,5,1,4,1,1,0,2,
       2,4,1,1,1,4,1,3,4,4,10,5,2,0,7,1,6,1,3,6,4,0,2,4,1,12,2,5,
       6,5,4,1,11,0,1,3,2,4,2,0,2,3,4,0,2,9,9,7,4,2,1,3,3,3,4,2,9,
       2,4,3,2,2,4,2,5,3,0,1,3,2,0,3,3,4,1,3,3,5,7,3,3,2,1,5,5,4,
       6,1,1,1,2,9,5,1,2,4,0,2,1,0,3,2,4,3,1,4,2,1,4,1,6,0,6,5,3,
       5,2,0,1,2,1,0,5,3,2,7,6,4,3,2,5,7,5,5,1,1,3,10,2,0,5,0,1,2,
       0,5,1,2,3,6,4,0,3,1,2,2,4,3,0,3,2,5,4,10,1,2,4,4,2,13,4,3,
       1,5,4,8,5,6,2,3,4,3,1,5,5,1,8,2,0,5,7,3,2,2,4,2,3,1,5,3,7,
       13,1,4,7,5,5,0,3,0,4,2,3,1,2,4,2,8,1,2,5,6,1,1,0,7,2,2,3,5,
       12,2,2,2,0,3,3,4,0,2,5,1,10,0,7,6,5,0,11,2,3,7,3,5,4,2,1,2,
       4,0,2,2,2,0,6,2,3,4,2,3,7,3,5,2,5,0,4,4,6,3,1,2,7,3,0,2,5,
       7,2,2,0,0,0,6,3,0,1,1,5,5,2,6,2,4,6,0,1,2,3,2,2,2,3,4,1,1,
       4,0,2,0,1,3,4,1,2,2,3,1,4,4,3,4,4,1,5,2,13,4,10,5,6,1,0,5,
       0,0,5,6,0,1,8,5,1,3,1,8,1,8,1,6,7,2,8,2,2,3,3,0,4,2,1,9,6,
       0,6,7,1,8,2,2,1,11,3,0,4,2,5,1,6,8,3,4,7,0,4,2,4,1,1,1,6,0,
       4,4,6,2,1,3,1,0,4,9,3,1,4,2,2,0,1};
   par = mle_nb(x5, 500);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->a-3.0057) < 1e-3);
   test_assert(par->pi == 0.0);
   test_assert(fabs(par->p-0.4854) < 1e-3);
   free(par);

   return;

}


void
test_mle_zinb
(void)
{

   // Cases checked by simulated annealing.
   
   zinb_par_t *par;
   
   // 0:53, 1:8, 2:9, 3:8, 4:4, 5:7, 6:3, 7:3, 8:1, 9:3, 10:1
   int x1[100] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,3,7,4,5,3,5,1,5,7,3,6,1,2,9,1,10,6,2,2,2,3,2,1,1,5,0,2,3,
      9,4,2,9,3,5,7,3,5,2,1,0,6,1,4,2,3,4,5,8,1 };

   par = mle_zinb(x1, 100);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->a - 3.6855) < 1e-3);
   test_assert(fabs(par->pi - 0.5110) < 1e-3);
   test_assert(fabs(par->p - 0.5044) < 1e-3);
   free(par);

   // 0:73, 1:7, 2:7, 3:3, 4:4, 5:2, 7:1, 8:2, 9:1
   int x2[100] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,
      1,0,0,8,1,1,3,0,0,8,2,4,1,2,0,3,2,0,0,4,0,0,3,1,0,2,0,0,5,
      7,0,0,2,4,0,2,1,0,0,0,0,0,0,0,1,2,9,0,4 };

   par = mle_zinb(x2, 100);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->a - 1.8251) < 1e-3);
   test_assert(fabs(par->pi - 0.3363) < 1e-3);
   test_assert(fabs(par->p - 0.4109) < 1e-3);
   free(par);

   return;

}

void
test_indexts
(void)
{
   int index[8];

   int ts_1[24] = {
      0, 0, 0,
      0, 0, 1,
      0, 1, 0,
      1, 0, 1,
      0, 0, 1,
      0, 0, 1,
      0, 1, 0,
      0, 0, 0,
   };

   int ts_2[27] = {
      -100, 0, 1,
      -1,   0, 100,
      -1,   0, 1,
      -1,   0, 100,
      -100, 0, 1,
      -1,   0, 1,
      -100, 0, 100,
       0,   0, 0,
      -100, 0, 100,
   };

   int expected_1[8] = {0,1,2,3,1,1,2,0};
   int expected_2[9] = {0,1,2,1,0,2,6,7,6};

   test_assert(indexts(8, 3, ts_1, index) == 0);
   for (int i = 0 ; i < 8 ; i++) {
      test_assert(index[i] == expected_1[i]);
   }

   test_assert(indexts(9, 3, ts_2, index) == 7);
   for (int i = 0 ; i < 8 ; i++) {
      test_assert(index[i] == expected_2[i]);
   }

}


void
test_fwdb
(void)
{
   int n = 3;
   int m = 2;
   const double Q[4] = {
      // transpose //
      0.8, 0.1,
      0.2, 0.9,
   };
   double init[2] = {.8, .2};
   double prob[6] = {
      // First series.
      0.1,      0.3,
      0.4,      0.8,
      // `fwd` can deal with emission probabilities in log space.
      log(0.2), log(0.7),
   };

   double *phi = malloc(m*n * sizeof(double));
   double *trans = malloc(m*m * sizeof(double));
   test_assert_critical(phi != NULL && trans != NULL);

   double ll = fwdb(m, n, Q, init, prob, phi, trans);

   double expected_loglik = -3.105547;
   double expected_alpha[6] = {
      0.5714286, 0.4285714,
      0.3333333, 0.6666667,
      0.1250000, 0.8750000,
   };
   double expected_phi[6] = {
      0.3571429, 0.6428571,
      0.1875000, 0.8125000,
      0.1250000, 0.8750000,
   };
   double expected_trans[4] = {
      // transpose //
      0.2714285, 0.0410714,
      0.2732143, 1.4142857,
   };

   test_assert(fabs(expected_loglik - ll) < 1e-6);
   for (int i = 0 ; i < m*n ; i++) {
      test_assert(fabs(expected_alpha[i] - prob[i]) < 1e-6);
      test_assert(fabs(expected_phi[i] - phi[i]) < 1e-6);
   }
   for (int i = 0 ; i < m*m ; i++) {
      test_assert(fabs(expected_trans[i] - trans[i]) < 1e-6);
   }

   free(phi);
   free(trans);

   return;

}


void
test_fwdb_NA
(void)
{

   int n = 3;
   int m = 2;
   const double Q[4] = {
      // transpose //
      0.8, 0.1,
      0.2, 0.9,
   };
   double init[2] = {.8, .2};
   double prob[6] = {
      // First series (NA in first position).
      log(-1.0),     0.3,
            0.4,     0.8,
      // `fwd` can deal log emission probabilities.
      log(0.2), log(0.7),
   };

   double *phi = malloc(m*n * sizeof(double));
   double *trans = malloc(m*m * sizeof(double));
   test_assert_critical(phi != NULL && trans != NULL);

   double ll = fwdb(m, n, Q, init, prob, phi, trans);

   double expected_loglik = -1.362578;
   double expected_alpha[12] = {
      0.8000000, 0.2000000,
      0.4925373, 0.5074627,
      0.1862500, 0.8137500,
   };
   double expected_phi[12] = {
      0.6250000, 0.3750000,
      0.3093750, 0.6906250,
      0.1862500, 0.8137500,
   };
   double expected_trans[4] = {
      // transpose //
      0.4650000, 0.0306250,
      0.4693750, 1.0350000,
   };

   test_assert(fabs(expected_loglik - ll) < 1e-6);
   for (int i = 0 ; i < m*n ; i++) {
      test_assert(fabs(expected_alpha[i] - prob[i]) < 1e-6);
      test_assert(fabs(expected_phi[i] - phi[i]) < 1e-6);
   }
   for (int i = 0 ; i < m*m ; i++) {
      test_assert(fabs(expected_trans[i] - trans[i]) < 1e-6);
   }

   free(phi);
   free(trans);

   return;

}


void
test_underflow
(void)
{
   int n = 3;
   int m = 2;
   const double Q[4] = {
      // transpose //
      0.8, 0.1,
      0.2, 0.9,
   };
   double init[2] = {.8, .2};
   double prob_1[6] = {
      // First series (underflow at first step).
      0.0,      0.0,
      0.4,      0.8,
      // `fwd` can deal log emission probabilities.
      log(0.2), log(0.7),
   };

   double *phi = malloc(m*n * sizeof(double));
   double *trans = malloc(m*m * sizeof(double));
   test_assert_critical(phi != NULL && trans != NULL);

   double ll_1 = fwdb(m, n, Q, init, prob_1, phi, trans);

   double expected_loglik_1 = -1.362578;
   double expected_alpha_1[6] = {
      0.8000000, 0.2000000,
      0.4925373, 0.5074627,
      0.1862500, 0.8137500,
   };
   double expected_phi_1[6] = {
      0.6250000, 0.3750000,
      0.3093750, 0.6906250,
      0.1862500, 0.8137500,
   };
   double expected_trans_1[4] = {
      // transpose //
      0.4650000, 0.0306250,
      0.4693750, 1.0350000,
   };

   test_assert(fabs(expected_loglik_1 - ll_1) < 1e-6);
   for (int i = 0 ; i < m*n ; i++) {
      test_assert(fabs(expected_alpha_1[i] - prob_1[i]) < 1e-6);
      test_assert(fabs(expected_phi_1[i] - phi[i]) < 1e-6);
   }
   for (int i = 0 ; i < m*m ; i++) {
      test_assert(fabs(expected_trans_1[i] - trans[i]) < 1e-6);
   }

   // Create the underflow at step 2.
   double prob_2[6] = {
      0.1,      0.3,
      0.0,      0.0,
      // `fwd` can deal log emission probabilities.
      log(0.2), log(0.7),
   };

   double ll_2 = fwdb(m, n, Q, init, prob_2, phi, trans);

   double expected_loglik_2 = -2.710553;
   double expected_alpha_2[6] = {
      0.5714286, 0.4285714,
      0.5000000, 0.5000000,
      0.1894737, 0.8105263,
   };
   double expected_phi_2[6] = {
      0.4451128, 0.5548872,
      0.3157895, 0.6842105,
      0.1894737, 0.8105263,
   };
   double expected_trans_2[4] = {
      // transpose //
      0.4571429, 0.0481203,
      0.3037594, 1.1909774,
   };

   test_assert(fabs(expected_loglik_2 - ll_2) < 1e-6);
   for (int i = 0 ; i < m*n ; i++) {
      test_assert(fabs(expected_alpha_2[i] - prob_2[i]) < 1e-6);
      test_assert(fabs(expected_phi_2[i] - phi[i]) < 1e-6);
   }
   for (int i = 0 ; i < m*m ; i++) {
      test_assert(fabs(expected_trans_2[i] - trans[i]) < 1e-6);
   }

   free(phi);
   free(trans);

   return;

}

void
test_block_fwdb
(void)
{

   size_t m = 2;
   size_t n = 6;
   double Q[4] = {
      // transpose //
      0.8, 0.1,
      0.2, 0.9,
   };
   double init[2] = {.8, .2};
   double prob[12] = {
      // First series.
      0.1,      0.3,
      0.4,      0.8,
      // `fwd` can deal log emission probabilities.
      log(0.2), log(0.7),
      // Second series (identical).
      0.1,      0.3,
      0.4,      0.8,
      0.2,      0.7,
   };

   double *phi = malloc(m*n * sizeof(double));
   double *trans = malloc(m*m * sizeof(double));
   test_assert_critical(phi != NULL && trans != NULL);

   unsigned int nblocks = 2;
   unsigned int size[2] = {3,3};

   double l = block_fwdb(m, nblocks, size, Q, init, prob, phi, trans);

   double expected_loglik = 2 * -3.105547;
   double expected_phi[12] = {
      0.3571429, 0.6428571,
      0.1875000, 0.8125000,
      0.1250000, 0.8750000,
      0.3571429, 0.6428571,
      0.1875000, 0.8125000,
      0.1250000, 0.8750000,
   };
   double expected_trans[4] = {
      // transpose //
      2 * 0.2714285, 2 * 0.0410714,
      2 * 0.2732143, 2 * 1.4142857,
   };

   test_assert(fabs(expected_loglik - l) < 1e-6);
   for (int i = 0 ; i < m*n ; i++) {
      test_assert(fabs(expected_phi[i] - phi[i]) < 1e-6);
   }
   for (int i = 0 ; i < m*m ; i++) {
      test_assert(fabs(expected_trans[i] - trans[i]) < 1e-6);
   }

   free(phi);
   free(trans);

   return;

}


void
test_block_fwdb_NA
(void)
{
   size_t m = 2;
   size_t n = 6;
   double Q[4] = {
      // transpose //
      0.8, 0.1,
      0.2, 0.9,
   };
   double init[2] = {.8, .2};
   double prob[12] = {
      // First series.
      log(-1.0), 0.3,
      0.4,       0.8,
      // `fwd` can deal log emission probabilities.
      log(0.2),  log(0.7),
      // Second series (identical).
      0.1,       0.3,
      log(-1.0), 0.8,
      0.2,       0.7,
   };

   double *phi = malloc(m*n * sizeof(double));
   double *trans = malloc(m*m * sizeof(double));
   test_assert_critical(phi != NULL && trans != NULL);

   unsigned int nblocks = 2;
   unsigned int size[2] = {3,3};

   double l = block_fwdb(m, nblocks, size, Q, init, prob, phi, trans);

   double expected_loglik = -1.362578 -2.710553;
   double expected_alpha[12] = {
      0.8000000, 0.2000000,
      0.4925373, 0.5074627,
      0.1862500, 0.8137500,
      0.5714286, 0.4285714,
      0.5000000, 0.5000000,
      0.1894737, 0.8105263,
   };
   double expected_phi[12] = {
      0.6250000, 0.3750000,
      0.3093750, 0.6906250,
      0.1862500, 0.8137500,
      0.4451128, 0.5548872,
      0.3157895, 0.6842105,
      0.1894737, 0.8105263,
   };
   double expected_trans[4] = {
      // transpose //
      0.4650000+0.4571429, 0.0306250+0.0481203,
      0.4693750+0.3037594, 1.0350000+1.1909774,
   };

   test_assert(fabs(expected_loglik - l) < 1e-6);
   for (int i = 0 ; i < m*n ; i++) {
      test_assert(fabs(expected_alpha[i] - prob[i]) < 1e-6);
      test_assert(fabs(expected_phi[i] - phi[i]) < 1e-6);
   }
   for (int i = 0 ; i < m*m ; i++) {
      test_assert(fabs(expected_trans[i] - trans[i]) < 1e-6);
   }

   free(phi);
   free(trans);

   return;

}


void
test_viterbi
(void)
{
   // -- Test 1, this is a test case I worked out manually with R -- //
   int m = 4;
   int n = 5;
   double Q_1[16] = {
      // transpose //
      log(0.80), log(0.05), log(0.05), log(0.05),
      log(0.10), log(0.90), log(0.30), log(0.35),
      log(0.05), log(0.05), log(0.50), log(0.45),
      log(0.05), log(0.00), log(0.25), log(0.15),
   };

   double init_1[4] = {log(0.05), log(0.10), log(0.05), log(0.80)};
   double prob_1[20] = {
      log(0.1), log(0.3), log(0.2), log(0.3),
      log(0.2), log(0.4), log(0.8), log(0.1),
      log(0.2), log(0.7), log(0.1), log(0.2),
      log(0.8), log(0.1), log(0.3), log(0.4),
      log(0.9), log(0.1), log(0.2), log(0.1),
   };

   int *path = malloc(n * sizeof(int));
   test_assert_critical(path != NULL);

   viterbi(
                    m,
                    n,
   (const double *) Q_1,
   (const double *) init_1,
   (const double *) prob_1,
                    path
   );

   //   state 0    state 1    state 2    state 3
   // -5.298317  -3.506558  -4.605170  -1.427116
   // -6.032287  -3.393229  -2.448768  -5.626821
   // -7.053938  -3.855265  -5.444500  -5.444500
   // -7.074140  -6.263210  -7.341620  -7.747085
   // -7.402645  -8.671156  -9.644205 -11.030499

   double expected_path_1[5] = {
      3,1,1,0,0,
   };

   for (int i = 0 ; i < n ; i++) {
      test_assert(expected_path_1[i] == path[i]);
   }

   free(path);
   path = NULL;

   //-- Test 2, this is a test using Wikipedia's Python example --//
   m = 2;
   n = 23;
   double Q_2[4] = {
      // transpose //
      log(0.7), log(0.3),
      log(0.4), log(0.6),
   };

   double init_2[2] = {log(0.6), log(0.4)};
   double prob_2[46] = {
      log(0.5), log(0.1), // normal
      log(0.4), log(0.3), // cold  
      log(0.1), log(0.6), // dizzy 
      log(0.1), log(0.6), // dizzy 
      log(0.5), log(0.1), // normal
      log(0.5), log(0.1), // normal
      log(0.5), log(0.1), // normal
      log(0.1), log(0.6), // dizzy 
      log(0.1), log(0.6), // dizzy 
      log(0.1), log(0.6), // dizzy 
      log(0.4), log(0.3), // cold  
      log(0.4), log(0.3), // cold  
      log(0.4), log(0.3), // cold  
      log(0.1), log(0.6), // dizzy 
      log(0.5), log(0.1), // normal
      log(0.5), log(0.1), // normal
      log(0.4), log(0.3), // cold  
      log(0.5), log(0.1), // normal
      log(0.4), log(0.3), // cold  
      log(0.1), log(0.6), // dizzy 
      log(0.1), log(0.6), // dizzy 
      log(0.5), log(0.1), // normal
      log(0.5), log(0.1), // normal
   };

   path = malloc(n * sizeof(int));
   test_assert_critical(path != NULL);

   viterbi(
                    m,
                    n,
   (const double *) Q_2,
   (const double *) init_2,
   (const double *) prob_2,
                    path
   );

   double expected_path_2[23] = {
      0,0,1,1,0,0,0,1,1,1,0,0,0,1,0,0,0,0,0,1,1,0,0,
   };

   for (int i = 0 ; i < n ; i++) {
      test_assert(expected_path_2[i] == path[i]);
   }

   free(path);

   return;

}


void
test_block_viterbi
(void)
{

   int m = 4;
   int n = 10;
   const double Q[16] = {
      // transpose //
      0.80, 0.05, 0.05, 0.05,
      0.10, 0.90, 0.30, 0.35,
      0.05, 0.05, 0.40, 0.45,
      0.05, 0.00, 0.25, 0.15,
   };
   double init[4] = {0.05, 0.10, 0.05, 0.80};
   double prob[40] = {
      // First series.
      0.1, 0.3, 0.2, 0.3,
      0.2, 0.4, 0.8, 0.1,
      0.2, 0.7, 0.1, 0.2,
      0.8, 0.1, 0.3, 0.4,
      0.9, 0.1, 0.2, 0.1,
      // Second series (identical).
      0.1, 0.3, 0.2, 0.3,
      0.2, 0.4, 0.8, 0.1,
      0.2, 0.7, 0.1, 0.2,
      0.8, 0.1, 0.3, 0.4,
      0.9, 0.1, 0.2, 0.1,
   };

   int *path = malloc(n * sizeof(int));
   test_assert_critical(path != NULL);

   int nblocks = 2;
   uint size[2] = {5,5};

   block_viterbi(
         m,
         nblocks,
         size,
         Q,
         init,
         prob,
         path
   );

   double expected_path[10] = {
      3,1,1,0,0,3,1,1,0,0,
   };

   for (int i = 0 ; i < n ; i++) {
      test_assert(expected_path[i] == path[i]);
   }


   const double log_Q[16] = {
      // transpose //
      log(0.80), log(0.05), log(0.05), log(0.05),
      log(0.10), log(0.90), log(0.30), log(0.35),
      log(0.05), log(0.05), log(0.40), log(0.45),
      log(0.05), log(0.00), log(0.25), log(0.15),
   };
   double log_init[4] = {log(0.05), log(0.10), log(0.05), log(0.80)};
   double log_prob[40] = {
      // First series.
      log(0.1), log(0.3), log(0.2), log(0.3),
      log(0.2), log(0.4), log(0.8), log(0.1),
      log(0.2), log(0.7), log(0.1), log(0.2),
      log(0.8), log(0.1), log(0.3), log(0.4),
      log(0.9), log(0.1), log(0.2), log(0.1),
      // Second series (identical).
      log(0.1), log(0.3), log(0.2), log(0.3),
      log(0.2), log(0.4), log(0.8), log(0.1),
      log(0.2), log(0.7), log(0.1), log(0.2),
      log(0.8), log(0.1), log(0.3), log(0.4),
      log(0.9), log(0.1), log(0.2), log(0.1),
   };

   block_viterbi(
         m,
         nblocks,
         size,
         log_Q,
         log_init,
         log_prob,
         path
   );

   for (int i = 0 ; i < n ; i++) {
      test_assert(expected_path[i] == path[i]);
   }

   free(path);

   return;

}


void
test_block_viterbi_NA
(void)
{

   int m = 4;
   int n = 10;
   const double Q[16] = {
      // transpose //
      0.80, 0.05, 0.05, 0.05,
      0.10, 0.90, 0.30, 0.35,
      0.05, 0.05, 0.40, 0.45,
      0.05, 0.00, 0.25, 0.15,
   };
   double init[4] = {0.05, 0.10, 0.05, 0.80};
   double prob[40] = {
      // First series.
     -1.0, 0.3, 0.2, 0.3,
      0.2, 0.4, 0.8, 0.1,
      0.2, 0.7, 0.1, 0.2,
      0.8, 0.1, 0.3, 0.4,
      0.9, 0.1, 0.2, 0.1,
      // Second series.
      0.1, 0.3, 0.2, 0.3,
     -1.0, 0.4, 0.8, 0.1,
      0.0, 0.0, 0.0, 0.0,
      0.8, 0.1, 0.3, 0.4,
      0.9, 0.1, 0.2, 0.1,
   };

   int *path = malloc(n * sizeof(int));
   test_assert_critical(path != NULL);

   int nblocks = 2;
   uint size[2] = {5,5};

   block_viterbi(
         m,
         nblocks,
         size,
         Q,
         init,
         prob,
         path
   );

   double expected_path[10] = {
      3,1,1,0,0,3,0,0,0,0,
   };

   for (int i = 0 ; i < n ; i++) {
      test_assert(expected_path[i] == path[i]);
   }

   // Set invalid initial probabilities.
   static double invalid_init_1[4] = {-1.0, 0.10, 0.05, 0.80};

   // Catch stderr.
   redirect_stderr();
   int exit_status = block_viterbi(
         m,
         nblocks,
         size,
         Q,
         invalid_init_1,
         prob,
         path
   );
   unredirect_stderr();

   test_assert(exit_status == -1);
   test_assert_stderr("mixed log/lin arguments in 'block_viterbi()'\n");

   static double invalid_init_2[4] = {.3, .3, .6, 0.0/0.0};

   // Catch stderr.
   redirect_stderr();
   exit_status = block_viterbi(
         m,
         nblocks,
         size,
         Q,
         invalid_init_2,
         prob,
         path
   );
   unredirect_stderr();

   test_assert(exit_status == -1);
   test_assert_stderr("invalid 'init' argument in 'block_viterbi()'\n");

   static double invalid_init_3[4] = {.3, .3, .6, .9};

   // Catch stderr.
   redirect_stderr();
   exit_status = block_viterbi(
         m,
         nblocks,
         size,
         Q,
         invalid_init_3,
         prob,
         path
   );
   unredirect_stderr();

   test_assert(exit_status == -1);
   test_assert_stderr("'init' is not a probability in 'block_viterbi()'\n");

   // Set invalid transition probabilities.
   static double invalid_Q_1[16] = {
      // transpose //
      -1.0, 0.05, 0.05, 0.05,
      0.10, 0.90, 0.30, 0.35,
      0.05, 0.05, 0.40, 0.45,
      0.05, 0.00, 0.25, 0.15,
   };

   // Catch stderr.
   redirect_stderr();
   exit_status = block_viterbi(
         m,
         nblocks,
         size,
         invalid_Q_1,
         init,
         prob,
         path
   );
   unredirect_stderr();

   test_assert(exit_status == -1);
   test_assert_stderr("mixed log/lin arguments in 'block_viterbi()'\n");

   static double invalid_Q_2[16] = {
      // transpose //
      0.80, 0.05, 0.05, 0.05,
      0.10, 0.90, 0.30, 0.35,
      0.05, 0.05, 0.40, 0.45,
      0.05, 0.00, 0.25, 0.0/0.0,
   };

   // Catch stderr.
   redirect_stderr();
   exit_status = block_viterbi(
         m,
         nblocks,
         size,
         invalid_Q_2,
         init,
         prob,
         path
   );
   unredirect_stderr();

   test_assert(exit_status == -1);
   test_assert_stderr("invalid 'Q' argument in 'block_viterbi()'\n");

   static double invalid_Q_3[16] = {
      // transpose //
      0.80, 0.05, 0.05, 0.05,
      0.10, 0.90, 0.30, 0.35,
      0.05, 0.05, 0.40, 0.45,
      0.05, 0.00, 0.25, 0.10,
   };

   // Catch stderr.
   redirect_stderr();
   exit_status = block_viterbi(
         m,
         nblocks,
         size,
         invalid_Q_3,
         init,
         prob,
         path
   );
   unredirect_stderr();

   test_assert(exit_status == -1);
   test_assert_stderr("'Q' is not stochastic in 'block_viterbi()'\n");

   free(path);

   return;

}


void
test_zinm_prob
(void)
{

   const uint m = 3;
   const uint n = 7;
   const uint r = 3;

   const double pi = 0.8;
   const double a = 1.2;
   // The constant C1 facilitates the definition of 'p'.
   const double C1 = 2.5;
   const double p[12] = {
      C1, 1.0, 1.0, 1.0, 
      C1, 1.0, 2.0, 1.0,
      C1, 1.0, 0.2, 0.1,
   };

   int y[21] = {
         1, 2, 2,
         0, 4, 2,
         1, 2, 2,
         1, 2, 0,
         0, 0, 0,
      // Underflow
      1500, 1, 2,
      // Negative values
        -1, 4, 2,
   };

   ChIP_t *ChIP = new_ChIP(m, 1, y, &n);

   double dummy[9] = {0};
   jahmm_t *jahmm = new_jahmm(r, ChIP);
   set_jahmm_par(jahmm, dummy, a, pi, p);

   int index[7] = {-1,-1,-1,-1,-1,-1,-1};
   double pem[21];
   indexts(n, r, y, index);

   //--               Test output type 0               --//
   int out = 0;
   redirect_stderr();
   zinm_prob(jahmm, index, out, pem);
   unredirect_stderr();

   double expected_pem_1[21] = {
      // Checked manually.
          7.713996e-05,      0.0001095280,      3.054426e-07,
          1.402544e-05,      6.740185e-05,      3.215185e-09,
          7.713996e-05,      0.0001095280,      3.054426e-07,
          0.0023334833,      0.0046275583,      0.0004410592,
          0.5105866396,      0.4541686425,      0.6840360201,
      -2563.1825314667,  -2813.7721384365,  -2013.2236637989,
          NAN,               NAN,               NAN,
   };

   for (int i = 0 ; i < 18 ; i++) {
      test_assert(fabs(expected_pem_1[i] - pem[i]) < 1e-9);
   }
   for (int i = 18 ; i < 21 ; i++) {
      test_assert(pem[i] != pem[i]);
   }

   // Test the warning message.
   test_assert_stderr("warning: renormalizing 'p'\n");

   //--               Test output type 1               --//
   out = 1;
   redirect_stderr();
   zinm_prob(jahmm, index, out, pem);
   unredirect_stderr();

   double expected_pem_2[21] = {
      // Checked manually.
      log(7.713996e-05), log(0.0001095280), log(3.054426e-07),
      log(1.402544e-05), log(6.740185e-05), log(3.215185e-09),
      log(7.713996e-05), log(0.0001095280), log(3.054426e-07),
      log(0.0023334833), log(0.0046275583), log(0.0004410592),
      log(0.5105866396), log(0.4541686425), log(0.6840360201),
      -2563.1825314667,  -2813.7721384365,  -2013.2236637989,
      NAN,          NAN,          NAN,
   };

   for (int i = 0 ; i < 18 ; i++) {
      test_assert(fabs(expected_pem_2[i] - pem[i]) < 1e-6);
   }
   for (int i = 18 ; i < 21 ; i++) {
      test_assert(pem[i] != pem[i]);
   }
   // Test the warning message.
   test_assert_stderr("warning: renormalizing 'p'\n");

   //--               Test output type 2               --//
   out = 2;
   redirect_stderr();
   zinm_prob(jahmm, index, out, pem);
   unredirect_stderr();

   double expected_pem_3[21] = {
      // Checked manually.
      7.713996e-05,      0.0001095280,      3.054426e-07,
      1.402544e-05,      6.740185e-05,      3.215185e-09,
      7.713996e-05,      0.0001095280,      3.054426e-07,
      0.0023334833,      0.0046275583,      0.0004410592,
      0.5105866396,      0.4541686425,      0.6840360201,
      0.0000000000,      0.0000000000,      0.0000000000,
      NAN,               NAN,               NAN,
   };

   for (int i = 0 ; i < 15 ; i++) {
      test_assert(fabs(expected_pem_3[i] - pem[i]) < 1e-8);
   }
   for (int i = 15 ; i < 18 ; i++) {
      test_assert(pem[i] == 0.0);
   }
   for (int i = 18 ; i < 21 ; i++) {
      test_assert(pem[i] != pem[i]);
   }
   // Test the warning message.
   test_assert_stderr("warning: renormalizing 'p'\n");

   free(ChIP);
   jahmm->ChIP = NULL;
   destroy_jahmm_all(jahmm);

   return;

}


void
test_bw_zinm
(void)
{
   int y[90] = {
      2,2,2,
      5,0,2,
      3,3,2,
      2,1,1,
      8,0,2,
      2,2,0,
      0,1,0,
      1,2,0,
      5,2,1,
      3,0,0,
      2,1,0,
      9,1,1,
      2,2,0,
      2,0,2,
      2,0,1,
      2,2,3,
      4,2,3,
      6,12,9,
      2,1,2,
      4,3,10,
      6,8,8,
      1,7,5,
      5,8,3,
      5,5,9,
      4,5,3,
      4,2,0,
      6,12,7,
      2,11,6,
      4,5,6,
      4,3,5,
   };

   unsigned size[2] = {15,15};
   ChIP_t *ChIP = new_ChIP(3, 2, y, size);

   double p[8] = {
      .3448276, .4137931, .1379310, .1034483,
      .1612903, .1935484, .3225806, .3225806,
   };
   double Q[4] = { .8, .2, .2, .8 };
   jahmm_t *jahmm = new_jahmm(2, ChIP);
   set_jahmm_par(jahmm, Q, 3.4, 1.0, p);
   bw_zinm(jahmm);

   double expected_Q[4] = {
      // transpose //
      1.000000, 0.000001,
      0.000000, 0.999999, 
   };
   for (size_t i = 0 ; i < 4 ; i++) {
      test_assert(fabs(jahmm->Q[i] - expected_Q[i]) < 1e-6);
   }

   free(ChIP);
   jahmm->ChIP = NULL;
   destroy_jahmm_all(jahmm);

   return;

}



void
test_update_trans
(void)
{
   double Q[9] = {0};
   double trans[9] = {
      // transpose //
      14, 56, 44,
      14, 57, 43,
      14, 58, 42,
   };

   update_trans(3, Q, trans);

   double expected_Q[9] = {
      // transpose //
      .3333, 0.3274, 0.3410,
      .3333, 0.3333, 0.3333,
      .3333, 0.3391, 0.3255,
   };

   for (size_t i = 0 ; i < 9 ; i++) {
      test_assert(fabs(Q[i] - expected_Q[i]) < 1e-3);
   }

   return;

}


void
test_read_file
(void)
{

   FILE *inputf = fopen("sample_file.txt", "r");
   test_assert_critical(inputf != NULL);
   ChIP_t *ChIP = read_file(inputf);
   fclose(inputf);

   test_assert_critical(ChIP != NULL);
   test_assert(ChIP->r == 3);

   test_assert(ChIP->nb == 1);
   test_assert(ChIP->size[0] == 5);

   int *y = ChIP->y;
   int expected_y[15] = {
      -1,-1,-1,
      -1,-1,-1,
      -1,-1,-1,
       0, 0, 0,
       2, 1, 0,
   };
   for (size_t i = 0 ; i < 15 ; i++) {
      test_assert(y[i] == expected_y[i]);
   }

   free(ChIP->y);
   free(ChIP);

   return;

}


void
test_eval_bw_f
(void)
//   double a,
//   double pi,
//   double p0,
//   double A,
//   double B,
//   double C,
//   double D,
//   double E
{
   return;
//   double term1 = (D + a*A) / p0;
//   double term2 = B * pi*a*pow(p0,a-1) / (pi*pow(p0,a)+1-pi);
//   return p0 + E/(term1 + term2) - 1.0 / C;
}

int
main(
   int argc,
   char **argv
)
{

   // Register test cases //
   const static test_case_t test_cases[] = {
      {"utils/new_histo", test_new_histo},
      {"utils/histo_push", test_histo_push},
      {"utils/compress_histo", test_compress_histo},
      {"utils/tabulate", test_tabulate},
      {"utils/indexts", test_indexts},
      {"ZINM/zinm_prob", test_zinm_prob},
      {"ZINM/eval_nb_f", test_eval_nb_f},
      {"ZINM/eval_nb_dfda", test_eval_nb_dfda},
      {"ZINM/eval_zinb_f", test_eval_zinb_f},
      {"ZINM/eval_zinb_g", test_eval_zinb_g},
      {"ZINM/eval_zinb_dfda", test_eval_zinb_dfda},
      {"ZINM/eval_zinb_dfdp", test_eval_zinb_dfdp},
      {"ZINM/eval_zinb_dgda", test_eval_zinb_dgda},
      {"ZINM/ll_zinb", test_ll_zinb},
      {"ZINM/nb_est_alpha", test_nb_est_alpha},
      {"ZINM/mle_nb", test_mle_nb},
      {"ZINM/mle_zinb", test_mle_zinb},
      {"HMM/fwdb", test_fwdb},
      {"HMM/fwdb (NAs)", test_fwdb_NA},
      {"HMM/fwdb (underflow)", test_underflow},
      {"HMM/block_fwdb", test_block_fwdb},
      {"HMM/block_fwdb (NAs)", test_block_fwdb_NA},
      {"HMM/viterbi", test_viterbi},
      {"HMM/block_viterbi", test_block_viterbi},
      {"HMM/block_viterbi (NAs)", test_block_viterbi_NA},
      {"BaumWelch/update_trans", test_update_trans},
      {"BaumWelch/bw_zinm", test_bw_zinm},
      {"JAHMM/eval_bw_f", test_eval_bw_f},
      {"JAHMM/read_file", test_read_file},
      {NULL, NULL}
   };

   return run_unittest(argc, argv, test_cases);

}
