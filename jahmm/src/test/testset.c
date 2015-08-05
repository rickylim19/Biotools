#include "test.h"

char error_buffer[1024];
int backup;


void
redirect_stderr_to
(char buffer[])
{
   // Flush stderr, redirect to /dev/null and set buffer.
   fflush(stderr);
   int temp = open("/dev/null", O_WRONLY);
   dup2(temp, STDERR_FILENO);
   memset(buffer, '\0', 1024 * sizeof(char));
   setvbuf(stderr, buffer, _IOFBF, 1024);
   close(temp);
   // Fill the buffer (needed for reset).
   fprintf(stderr, "fill the buffer");
   fflush(stderr);
}

void
unredirect_sderr
(void)
{ 
   fflush(stderr);
   dup2(backup, STDERR_FILENO);
   setvbuf(stderr, NULL, _IONBF, 0);
}


// ----------------------------  utils.c ---------------------------- //

void
test_mean
(void)
{

   // -- First test case -- //
   int yz_1[8] = {
      1, 2,
      3, 4,
      5, 6,
      7, 8,
   };

   g_assert_cmpfloat(mean(yz_1, 4, 2), ==, 4.0);
   g_assert_cmpfloat(mean(yz_1+1, 4, 2), ==, 5.0);

   // -- Second test case -- //
   int NA = (int) log(-1);
   int yz_2[8] = {
       1, NA,
       3, NA,
      NA, NA,
       8, NA,
   };

   g_assert_cmpfloat(mean(yz_2, 4, 2), ==, 4.0);
   g_assert_cmpfloat(mean(yz_2+1, 4, 2), !=, mean(yz_2+1, 4, 2));

}

void
test_histsum
(void)
{

   // -- First test case -- //
   int yz_1[8] = {
      1, 3,
      2, 2,
      3, 1,
      4, 0,
   };

   int expected_hist_1[5] = {0,0,0,0,4};
   int *hist_1 = histsum(yz_1, 4, 2);

   for (int i = 0 ; i < 5 ; i++) {
      g_assert_cmpint(hist_1[i], ==, expected_hist_1[i]);
   }

   free(hist_1);

   // -- Second test case -- //
   int yz_2[8] = {
       1, -3,
       2,  2,
      -3, -1,
       4,  0,
   };

   int expected_hist_2[5] = {0,0,0,0,2};
   int *hist_2 = histsum(yz_2, 4, 2);

   for (int i = 0 ; i < 5 ; i++) {
      g_assert_cmpint(hist_2[i], ==, expected_hist_2[i]);
   }

   free(hist_2);

   // -- Third test case -- //
   int yz_3[8] = {
       1, -3,
       2, -2,
      -3, -1,
       0,  1024,
   };
   int *hist_3 = histsum(yz_3, 4, 2);

   g_assert_cmpint(hist_3[1024], ==, 1);
   for (int i = 0 ; i < 1024 ; i++) {
      g_assert_cmpint(hist_3[i], ==, 0);
   }

   free(hist_3);

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

   int ts_2[24] = {
      -100, 0, 1,
      -1,   0, 100,
      -1,   0, 1,
      -1,   0, 100,
      -100, 0, 1,
      -1,   0, 1,
      -100, 0, 100,
      -100, 0, 100,
   };

   int expected_1[8] = {0,1,2,3,1,1,2,0};
   int expected_2[8] = {0,1,2,1,0,2,6,6};

   indexts(8, 3, ts_1, index);
   for (int i = 0 ; i < 8 ; i++) {
      g_assert_cmpint(index[i], ==, expected_1[i]);
   }

   indexts(8, 3, ts_2, index);
   for (int i = 0 ; i < 8 ; i++) {
      g_assert_cmpint(index[i], ==, expected_2[i]);
   }

}


// -------------------------  mnmultinom.c ------------------------- //

void
test_mnmultinom_prob
(void)
{

   int m = 3;
   int n = 6;
   int r = 3;

   double t = 0.8;
   double a = 1.2;
   // 'C1' and 'C2' simplify the definition of 'p' and 'q'.
   double C1 = 0.7692308;
   double C2 = 2.5000000;
   double p[12] = {
      C1/(C1+3.0), 1.0/(C1+3.0), 1.0/(C1+3.0), 1.0/(C1+3.0), 
      C1/(C1+4.0), 1.0/(C1+4.0), 2.0/(C1+4.0), 1.0/(C1+4.0),
      C1/(C1+1.3), 1.0/(C1+1.3), 0.2/(C1+1.3), 0.1/(C1+1.3),
   };
   double q[12] = {
      C2/(C2+3.0),   1.0/(C2+3.0),   1.0/(C2+3.0),   1.0/(C2+3.0), 
      C2/(C2+4.0),   1.0/(C2+4.0),   2.0/(C2+4.0),   1.0/(C2+4.0),
      // This line is not properly normalized. It should
      // trigger a warning but not cause failure.
      C2/(C2+1.3)*2, 1.0/(C2+1.3)*2, 0.2/(C2+1.3)*2, 0.1/(C2+1.3)*2,
   };

   int yz[18] = {
         1, 2, 2,
         0, 4, 2,
         1, 2, 2,
         1, 2, 0,
         // Underflow (should give 0.0).
      1500, 1, 2,
         // Negative values (should give NA).
        -1, 4, 2,
   };

   int index[6] = {-1,-1,-1,-1,-1,-1};
   double pem[18];
   indexts(n, r, (const int *) yz, index);

   //--               Test output type 0               --//
   int out = 0;
   redirect_stderr_to(error_buffer);
   mnmultinom_prob(&m, &n, &r, yz, &t, &a, p, q, index, &out, pem);
   unredirect_sderr();

   double expected_pem_1[18] = {
      // Checked manually.
          0.0001715956,      0.0001671296,      2.633868e-06,
          4.423732e-05,      0.0001352811,      5.037707e-08,
          0.0001715956,      0.0001671296,      2.633868e-06,
          0.0026853860,      0.0042287120,      0.0011898030,
      -1996.41726182421, -2349.68356163581, -1100.57214790818,
      NAN,          NAN,          NAN,
   };
   char expected_warning[] = "warning: renormalizing 'p' and/or 'q'\n";

   for (int i = 0 ; i < 15 ; i++) {
      g_assert_cmpfloat(fabs(expected_pem_1[i] - pem[i]), <, 1e-8);
   }
   for (int i = 15 ; i < 18 ; i++) {
      g_assert(pem[i] != pem[i]);
   }
   // Test the warning message.
   g_assert_cmpstr(error_buffer, ==, expected_warning);

   //--               Test output type 1               --//
   out = 1;
   redirect_stderr_to(error_buffer);
   mnmultinom_prob(&m, &n, &r, yz, &t, &a, p, q, index, &out, pem);
   unredirect_sderr();

   double expected_pem_2[18] = {
      // Checked manually.
      log(0.0001715956), log(0.0001671296), log(2.633868e-06),
      log(4.423732e-05), log(0.0001352811), log(5.037707e-08),
      log(0.0001715956), log(0.0001671296), log(2.633868e-06),
      log(0.002685386),  log(0.004228712),  log(0.001189803),
      -1996.41726182421, -2349.68356163581, -1100.57214790818,
      NAN,          NAN,          NAN,
   };

   for (int i = 0 ; i < 15 ; i++) {
      g_assert_cmpfloat(fabs(expected_pem_2[i] - pem[i]), <, 1e-6);
   }
   for (int i = 15 ; i < 18 ; i++) {
      g_assert(pem[i] != pem[i]);
   }
   // Test the warning message.
   g_assert_cmpstr(error_buffer, ==, expected_warning);

   //--               Test output type 2               --//
   out = 2;
   redirect_stderr_to(error_buffer);
   mnmultinom_prob(&m, &n, &r, yz, &t, &a, p, q, index, &out, pem);
   unredirect_sderr();

   double expected_ratio[18] = {
      // Checked manually.
      0.9100910, 0.8689304, 0.9768065,
      0.9365900, 0.9003529, 0.9872355,
      0.9100910, 0.8689304, 0.9768065,
      0.8262087, 0.7811363, 0.9258601,
      1.0000000, 1.0000000, 1.0000000,
      NAN,       NAN,       NAN,
   };

   for (int i = 0 ; i < 15 ; i++) {
      g_assert_cmpfloat(fabs(expected_ratio[i] - pem[i]), <, 1e-7);
   }
   for (int i = 15 ; i < 18 ; i++) {
      g_assert(pem[i] != pem[i]);
   }
   // Test the warning message.
   g_assert_cmpstr(error_buffer, ==, expected_warning);

   //--               Test output type 3               --//
   out = 3;
   redirect_stderr_to(error_buffer);
   mnmultinom_prob(&m, &n, &r, yz, &t, &a, p, q, index, &out, pem);
   unredirect_sderr();

   double expected_pem_3[18] = {
      // Checked manually.
      0.0001715956, 0.0001671296, 2.633868e-06,
      4.423732e-05, 0.0001352811, 5.037707e-08,
      0.0001715956, 0.0001671296, 2.633868e-06,
      0.002685386,  0.004228712,  0.001189803,
      0.0000000000, 0.0000000000, 0.0000000000,
      NAN,          NAN,          NAN,
   };

   for (int i = 0 ; i < 12 ; i++) {
      g_assert_cmpfloat(fabs(expected_pem_3[i] - pem[i]), <, 1e-8);
   }
   for (int i = 12 ; i < 15 ; i++) {
      g_assert_cmpfloat(pem[i], ==, 0.0);
   }
   for (int i = 15 ; i < 18 ; i++) {
      g_assert(pem[i] != pem[i]);
   }
   // Test the warning message.
   g_assert_cmpstr(error_buffer, ==, expected_warning);


   //--         Test output type 1 with constant terms         --//
   out = 1 + 8;
   redirect_stderr_to(error_buffer);
   mnmultinom_prob(&m, &n, &r, yz, &t, &a, p, q, index, &out, pem);
   unredirect_sderr();

   double expected_pem_4[18] = {
      // Checked manually.
      log(0.0001715956) + 3.83137851672, log(0.0001671296) + 3.83137851672,
            log(2.633868e-06) + 3.83137851672,
      log(4.423732e-05) + 3.17102115898, log(0.0001352811) + 3.17102115898,
            log(5.037707e-08) + 3.17102115898,
      log(0.0001715956) + 3.83137851672, log(0.0001671296) + 3.83137851672,
            log(2.633868e-06) + 3.83137851672,
      log(0.002685386) + 1.4407825464,  log(0.004228712) + 1.4407825464,
            log(0.001189803) + 1.4407825464,
      -1996.41726182421 + 22.799008469, -2349.68356163581 + 22.799008469,
            -1100.57214790818 + 22.799008469,
      NAN,          NAN,          NAN,
   };

   for (int i = 0 ; i < 15 ; i++) {
      g_assert_cmpfloat(fabs(expected_pem_4[i] - pem[i]), <, 1e-6);
   }
   for (int i = 15 ; i < 18 ; i++) {
      g_assert(pem[i] != pem[i]);
   }
   // Test the warning message.
   g_assert_cmpstr(error_buffer, ==, expected_warning);

   //--         Test output type 3 with constant terms         --//
   out = 3 + 8;
   redirect_stderr_to(error_buffer);
   mnmultinom_prob(&m, &n, &r, yz, &t, &a, p, q, index, &out, pem);
   unredirect_sderr();

   double expected_pem_5[18] = {
      // Checked manually.
      0.0001715956 * 46.12608,  0.0001671296 * 46.12608,
            2.633868e-06 * 46.12608,
      4.423732e-05 * 23.831808, 0.0001352811 * 23.831808,
            5.037707e-08 * 23.831808, 
      0.0001715956 * 46.12608,  0.0001671296 * 46.12608,
            2.633868e-06 * 46.12608,
      0.002685386 * 4.224,      0.004228712 * 4.224,
            0.001189803 * 4.224,
      0.000000,     0.000000,     0.000000,
      NAN,          NAN,          NAN,
   };

   for (int i = 0 ; i < 12 ; i++) {
      g_assert_cmpfloat(fabs(expected_pem_5[i] - pem[i]), <, 1e-6);
   }
   for (int i = 12 ; i < 15 ; i++) {
      g_assert_cmpfloat(pem[i], ==, 0);
   }
   for (int i = 15 ; i < 18 ; i++) {
      g_assert(pem[i] != pem[i]);
   }
   // Test the warning message.
   g_assert_cmpstr(error_buffer, ==, expected_warning);

   return;

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

   g_assert_cmpfloat(fabs(expected_loglik - ll), <, 1e-6);
   for (int i = 0 ; i < m*n ; i++) {
      g_assert_cmpfloat(fabs(expected_alpha[i] - prob[i]), <, 1e-6);
      g_assert_cmpfloat(fabs(expected_phi[i] - phi[i]), <, 1e-6);
   }
   for (int i = 0 ; i < m*m ; i++) {
      g_assert_cmpfloat(fabs(expected_trans[i] - trans[i]), <, 1e-6);
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

   // Make sure we catch stderr.
   redirect_stderr_to(error_buffer);

   double ll = fwdb(m, n, Q, init, prob, phi, trans);

   unredirect_sderr();

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

   g_assert_cmpfloat(fabs(expected_loglik - ll), <, 1e-6);
   for (int i = 0 ; i < m*n ; i++) {
      g_assert_cmpfloat(fabs(expected_alpha[i] - prob[i]), <, 1e-6);
      g_assert_cmpfloat(fabs(expected_phi[i] - phi[i]), <, 1e-6);
   }
   for (int i = 0 ; i < m*m ; i++) {
      g_assert_cmpfloat(fabs(expected_trans[i] - trans[i]), <, 1e-6);
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

   g_assert_cmpfloat(fabs(expected_loglik_1 - ll_1), <, 1e-6);
   for (int i = 0 ; i < m*n ; i++) {
      g_assert_cmpfloat(fabs(expected_alpha_1[i] - prob_1[i]), <, 1e-6);
      g_assert_cmpfloat(fabs(expected_phi_1[i] - phi[i]), <, 1e-6);
   }
   for (int i = 0 ; i < m*m ; i++) {
      g_assert_cmpfloat(fabs(expected_trans_1[i] - trans[i]), <, 1e-6);
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

   g_assert_cmpfloat(fabs(expected_loglik_2 - ll_2), <, 1e-6);
   for (int i = 0 ; i < m*n ; i++) {
      g_assert_cmpfloat(fabs(expected_alpha_2[i] - prob_2[i]), <, 1e-6);
      g_assert_cmpfloat(fabs(expected_phi_2[i] - phi[i]), <, 1e-6);
   }
   for (int i = 0 ; i < m*m ; i++) {
      g_assert_cmpfloat(fabs(expected_trans_2[i] - trans[i]), <, 1e-6);
   }

   free(phi);
   free(trans);

   return;

}

void
test_block_fwdb
(void)
{

   int m = 2;
   int n = 6;
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

   int nblocks = 2;
   int size[2] = {3,3};

   double ll;
   block_fwdb(&m, &nblocks, size, Q, init, prob, phi, trans, &ll);

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

   g_assert_cmpfloat(fabs(expected_loglik - ll), <, 1e-6);
   for (int i = 0 ; i < m*n ; i++) {
      g_assert_cmpfloat(fabs(expected_phi[i] - phi[i]), <, 1e-6);
   }
   for (int i = 0 ; i < m*m ; i++) {
      g_assert_cmpfloat(fabs(expected_trans[i] - trans[i]), <, 1e-6);
   }

   free(phi);
   free(trans);

   return;

}


void
test_block_fwdb_NA
(void)
{
   int m = 2;
   int n = 6;
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

   int nblocks = 2;
   int size[2] = {3,3};

   double ll;
   block_fwdb(&m, &nblocks, size, Q, init, prob, phi, trans, &ll);

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

   g_assert_cmpfloat(fabs(expected_loglik - ll), <, 1e-6);
   for (int i = 0 ; i < m*n ; i++) {
      g_assert_cmpfloat(fabs(expected_alpha[i] - prob[i]), <, 1e-6);
      g_assert_cmpfloat(fabs(expected_phi[i] - phi[i]), <, 1e-6);
   }
   for (int i = 0 ; i < m*m ; i++) {
      g_assert_cmpfloat(fabs(expected_trans[i] - trans[i]), <, 1e-6);
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
      g_assert_cmpint(expected_path_1[i], ==, path[i]);
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
      g_assert_cmpint(expected_path_2[i], ==, path[i]);
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
      0.05, 0.05, 0.50, 0.45,
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
   int nblocks = 2;
   static int size[2] = {5,5};
   int nolog = 0;

   block_viterbi(
         &m,
         &nblocks,
         size,
         Q,
         init,
         prob,
         &nolog,
         path
   );

   double expected_path[10] = {
      3,1,1,0,0,3,1,1,0,0,
   };

   for (int i = 0 ; i < n ; i++) {
      g_assert(expected_path[i] == path[i]);
   }


   int yeslog = 1;
   const double log_Q[16] = {
      // transpose //
      log(0.80), log(0.05), log(0.05), log(0.05),
      log(0.10), log(0.90), log(0.30), log(0.35),
      log(0.05), log(0.05), log(0.50), log(0.45),
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
         &m,
         &nblocks,
         size,
         log_Q,
         log_init,
         log_prob,
         &yeslog,
         path
   );

   for (int i = 0 ; i < n ; i++) {
      g_assert_cmpint(expected_path[i], ==, path[i]);
   }

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
      0.05, 0.05, 0.50, 0.45,
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
   int nblocks = 2;
   static int size[2] = {5,5};
   int nolog = 0;

   block_viterbi(
         &m,
         &nblocks,
         size,
         Q,
         init,
         prob,
         &nolog,
         path
   );

   double expected_path[10] = {
      3,1,1,0,0,3,0,0,0,0,
   };

   for (int i = 0 ; i < n ; i++) {
      g_assert_cmpint(expected_path[i], ==, path[i]);
   }

   // Set invalid initial probabilities.
   static double invalid_init[4] = {-1.0, 0.10, 0.05, 0.80};

   // Catch stderr.
   redirect_stderr_to(error_buffer);

   int exit_status = block_viterbi(
         &m,
         &nblocks,
         size,
         Q,
         invalid_init,
         prob,
         &nolog,
         path
   );

   unredirect_sderr();

   char expected_warning_2[] =
      "invalid 'init' parameter in 'block_viterbi'\n";

   g_assert_cmpint(exit_status, ==, -1);

   // Check the error message.
   g_assert_cmpstr(error_buffer, ==, expected_warning_2);

   // Set invalid transition probabilities.
   static double invalid_Q[16] = {
      // transpose //
      -1.0, 0.05, 0.05, 0.05,
      0.10, 0.90, 0.30, 0.35,
      0.05, 0.05, 0.50, 0.45,
      0.05, 0.00, 0.25, 0.15,
   };

   // Catch stderr.
   redirect_stderr_to(error_buffer);

   exit_status = block_viterbi(
         &m,
         &nblocks,
         size,
         invalid_Q,
         init,
         prob,
         &nolog,
         path
   );

   unredirect_sderr();

   char expected_warning_3[] =
      "invalid 'Q' parameter in 'block_viterbi'\n";

   g_assert_cmpint(exit_status, ==, -1);

   // Check the error message.
   g_assert_cmpstr(error_buffer, ==, expected_warning_3);

   free(path);

   return;

}


void
test_perf_mnmultinom_prob
(void)
{
   // -- Open input file (the code is hopelessly not portable) -- //
   FILE *f = fopen("RXRA.txt", "r");

   int n = 1031884;
   int i = 0;

   char seqname[32];
   int *yz = malloc (3*n * sizeof(int));
   memset(yz, (int) -1, 3*n * sizeof(int));

   // `getline` is available because we define _GNU_SOURCE.
   char *line = NULL;
   size_t len = 0;
   ssize_t read;
   // Discard header (the 'if' turns off unused variable warnings).
   if (getline(&line, &len, f));
   while ((read = getline(&line, &len, f)) != -1) {
      sscanf(line, "%s\t%d\t%d\t%d\n", seqname, 
         yz+i, yz+i+1, yz+i+2);
      i += 3;
   }
   free(line);
   fclose(f);

   // --               Run the performance test               -- //
   int m = 3;
   int r = 3;

   double t = 0.91;
   double a = 2.844083;

   double p[12] = {
      // transpose //
      0.05155538, 0.75760266, 0.11489932, 0.07594263,
      0.05066482, 0.74451597, 0.15030383, 0.05451538,
      0.04539102, 0.66701779, 0.15180765, 0.13578353,
   };
   double q[12] = {
      // transpose //
      0.54836008, 0.38167944, 0.04426865, 0.02569183,
      0.52542413, 0.36571515, 0.07734838, 0.03151234,
      0.47888700, 0.33332350, 0.16109890, 0.02669060,
   };
   double Q[9] = {
      // transpose //
      0.8712810672, 0.0004510562, 0.4013737724,
      0.0016887750, 0.9669194640, 0.3386870900,
      0.1270301600, 0.0326294800, 0.2599391400,
   };
   double init[3] = {.33333, .33333, .33333};
   double trans[9];

   int *index = malloc(n * sizeof(int));
   double *pem = malloc(n*m * sizeof(double));
   double *phi = malloc(n*m * sizeof(double));
   indexts(n, r, (const int *) yz, index);

   int out = 0;
   double ll;
   int nblocks = 24;
   int size[24] = {
      83083, 45178, 45002, 44617, 38389, 35783, 34177,
      30118, 27065, 26025, 19709, 81066, 21008, 16043,
      17101, 66007, 63718, 60305, 57038, 53046, 48788,
      47071, 51756, 19791,
   };
   redirect_stderr_to(error_buffer);
   mnmultinom_prob(&m, &n, &r, yz, &t, &a, p, q, index, &out, pem);
   block_fwdb(&m, &nblocks, size, Q, init, pem, phi, trans, &ll);
   unredirect_sderr();

   free(yz);
   free(index);
   free(pem);
   free(phi);

}

void
test_params
(void)
{

   const int m = 3;
   const int r = 3;

   double Q[9] = {
      0.970078, 0.007204, 0.025482,
      0.010285, 0.961420, 0.017606,
      0.019637, 0.031376, 0.956912,
   };
   double p[12] = {
      .0548, .8062, .1304, .0086,
      .0459, .6758, .1569, .1214,
      .0399, .5865, .2199, .1537,
   };
   double q[12] = {
      .5331, .3711, .0611, .0347,
      .5387, .3750, .0245, .0618,
      .4807, .3345, .0933, .0915,
   };

   // Test `params_new`.
   params *par = params_new(m, r);
   params *try = params_new(m, r);
   g_assert(par != NULL);
   g_assert(try != NULL);
   g_assert_cmpint(par->m, ==, m);
   g_assert_cmpint(par->r, ==, r);
   g_assert_cmpint(try->m, ==, m);
   g_assert_cmpint(try->r, ==, r);

   // Test `params_set`.
   params_set(par, 0.9172, 2.8440, Q, p, q);
   g_assert_cmpfloat(par->t, ==, 0.9172);
   g_assert_cmpfloat(par->a, ==, 2.8440);
   for (int i = 0 ; i < m*m ; i++) g_assert(par->Q[i] == Q[i]);
   for (int i = 0 ; i < m*(r+1) ; i++) {
      g_assert_cmpfloat(par->p[i], ==, p[i]);
      g_assert_cmpfloat(par->q[i], ==, q[i]);
   }

   // Test `params_cpy`.
   params_cpy(try, par);
   g_assert_cmpfloat(try->t, ==, 0.9172);
   g_assert_cmpfloat(try->a, ==, 2.8440);
   for (int i = 0 ; i < m*m ; i++) {
      g_assert_cmpfloat(try->Q[i], ==, Q[i]);
   }
   for (int i = 0 ; i < m*(r+1) ; i++) {
      g_assert_cmpfloat(try->p[i], ==, p[i]);
      g_assert_cmpfloat(try->q[i], ==, q[i]);
   }
   
   // Test `params_change`.
   for (int i = 0 ; i < 262144 ; i++) {
      params_change(try, par, i % 3);
      for (int j = 0 ; j < m ; j++) {
         double sumQ = 0.0;
         for (int k = 0 ; k < m ; k++) {
            sumQ += try->Q[j+k*m];
            g_assert_cmpfloat(try->Q[j+k*m], <, 1);
            g_assert_cmpfloat(try->Q[j+k*m], >, 0);
         }
         g_assert_cmpfloat(fabs(sumQ - 1.0), <, 1e-6);
         double sump = 0.0;
         double sumq = 0.0;
         for (int k = 0 ; k < r+1 ; k++) {
            sump += try->p[k+j*(r+1)];
            sumq += try->q[k+j*(r+1)];
            g_assert_cmpfloat(try->p[k+j*(r+1)], <, 1);
            g_assert_cmpfloat(try->p[k+j*(r+1)], >, 0);
            g_assert_cmpfloat(try->q[k+j*(r+1)], <, 1);
            g_assert_cmpfloat(try->q[k+j*(r+1)], >, 0);
         }
         g_assert_cmpfloat(fabs(sump - 1.0), <, 1e-6);
         g_assert_cmpfloat(fabs(sumq - 1.0), <, 1e-6);
      }
   }

   params_destroy(try);
   params_destroy(par);

   return;

}

void
test_loglik
(void)
{

   int m = 3;
   int n = 6;
   int r = 3;

   double t = 0.8;
   double a = 1.2;

   double Q[9] = {
      // transpose //
      0.7, 0.1, 0.9,
      0.2, 0.5, 0.1,
      0.1, 0.4, 0.0,
   };

   // 'C1' and 'C2' simplify the definition of 'p' and 'q'.
   double C1 = 0.7692308;
   double C2 = 2.5000000;
   double p[12] = {
      C1/(C1+3.0), 1.0/(C1+3.0), 1.0/(C1+3.0), 1.0/(C1+3.0), 
      C1/(C1+4.0), 1.0/(C1+4.0), 2.0/(C1+4.0), 1.0/(C1+4.0),
      C1/(C1+1.3), 1.0/(C1+1.3), 0.2/(C1+1.3), 0.1/(C1+1.3),
   };
   double q[12] = {
      C2/(C2+3.0),   1.0/(C2+3.0),   1.0/(C2+3.0),   1.0/(C2+3.0), 
      C2/(C2+4.0),   1.0/(C2+4.0),   2.0/(C2+4.0),   1.0/(C2+4.0),
      // This line is not properly normalized. It should
      // trigger a warning but not cause failure.
      C2/(C2+1.3)*2, 1.0/(C2+1.3)*2, 0.2/(C2+1.3)*2, 0.1/(C2+1.3)*2,
   };

   params *par = params_new(m, r);
   params_set(par, t, a, Q, p, q);

   int yz[18] = {
         1, 2, 2,
         0, 4, 2,
         1, 2, 2,
         1, 2, 0,
         // Underflow (should give 0.0).
      1500, 1, 2,
         // Negative values (should give NA).
        -1, 4, 2,
   };

   int index[6] = {-1,-1,-1,-1,-1,-1};
   double pem[18];
   indexts(n, r, (const int *) yz, index);

   redirect_stderr_to(error_buffer);
   double ll = loglik(n, par, yz, index, pem);
   unredirect_sderr();

   double expected_ll = -1135.57472829;
   double expected_alpha[18] = {
      // Checked manually.
          0.502683584781,  0.489600586793,  0.007715828424,
          0.278067313112,  0.721741573153,  0.000191113733,
          0.394077474900,  0.598752057680,  0.007170467419,
          0.322075687622,  0.561610874783,  0.116313437595,
          0.000000000000,  0.000000000000,  1.000000000000,
          0.900000000000,  0.100000000000,  0.000000000000,
   };   

   char expected_warning[] = "fill the buffer";

   g_assert_cmpfloat(abs(ll - expected_ll), <, 1e-6);

   for (int i = 0 ; i < 18 ; i++) {
      g_assert_cmpfloat(fabs(expected_alpha[i] - pem[i]), <, 1e-6);
   }

   // Test the warning message.
   g_assert_cmpstr(error_buffer, ==, expected_warning);

   params_destroy(par);

   return;
}


int
main(
   int argc,
   char **argv
)
{

   // Store 'stderr', file descriptor.
   backup = dup(STDERR_FILENO);

   g_test_init(&argc, &argv, NULL);
   // utils.c //
   g_test_add_func("/utils/indexts", test_indexts);
   g_test_add_func("/utils/histsum", test_histsum);
   g_test_add_func("/utils/mean", test_mean);
   // mnmultinom.c //
   g_test_add_func("/mnmultinom_prob", test_mnmultinom_prob);
   // hmm.c //
   g_test_add_func("/fwd-bwd/fwdb", test_fwdb);
   g_test_add_func("/fwd-bwd/fwdb_NA", test_fwdb_NA);
   g_test_add_func("/fwd-bwd/underflow", test_underflow);
   g_test_add_func("/fwd-bwd/block_fwdb_NA", test_block_fwdb);
   g_test_add_func("/fwd-bwd/block_fwdb_NA", test_block_fwdb_NA);
   g_test_add_func("/viterbi/viterbi", test_viterbi);
   g_test_add_func("/viterbi/block_viterbi", test_block_viterbi);
   g_test_add_func("/viterbi/block_viterbi_NA", test_block_viterbi_NA);
   // pso.c //
   g_test_add_func("/pso/params_change", test_params);
   g_test_add_func("/pso/loglik", test_loglik);

   // performance profiling //
   if (g_test_perf()) {
      g_test_add_func("/prof/e_step", test_perf_mnmultinom_prob);
   }

   int g_test_result = g_test_run();
   close(backup);

   return g_test_result;

}
