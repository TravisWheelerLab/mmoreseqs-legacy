/*******************************************************************************
 *  @file misc.h
 *  @brief Miscellaneous Helper Functions.
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/

/* external imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/* local imports (after struct declarations) */
#include "structs.h"
#include "misc.h"
#include "hmm_parser.h"

/* Global Vars */
/* table of logsum values */
static float logsum_lookup[LOGSUM_TBL];
bool logsum_initialized = false;

/* TODO: inline small functions? */

/* Get max value of two floats */
// static inline
float calc_Max (float x, float y)
{
   if (x > y)
   {
      return x;
   }
   return y;
}

/* Get min value of two floats */
// static inline
float calc_Min (float x, float y)
{
   if (x < y)
   {
      return x;
   }
   return y;
}

/* initialize the logsum table */
void init_Logsum ()
{
   if (logsum_initialized == false)
   {
      logsum_initialized = true;
      int i;
      for (i = 0; i < LOGSUM_TBL; i++)
      {
         logsum_lookup[i] = log(1. + exp((double) - i / LOGSUM_SCALE));
      }
   }

   // /* DAVID RICH EDIT */
   // printf("LOGSUM TABLE:\n");
   // FILE *fp = fopen("logsum.mine.txt", "w+");
   // for (int i = 0; i < LOGSUM_TBL; i++)
   // {
   //    fprintf(fp, "%d -> %f\n", i, logsum_lookup[i]);
   // }
   // fclose(fp);
   // exit(0);

}

/* Takes two logscale numbers and returns the log of their real sum (approx) */
/* ans = log( exp(x) + exp(y) ) */
// static inline
float calc_Logsum (float x, float y)
{
   float max, min;

   if (x > y)
   {
      max = x;
      min = y;
   }
   else
   {
      max = y;
      min = x;
   }

   return (min == -INF || (max - min) >= 15.7f) ?
          max : max + logsum_lookup[ (int)((max - min) * LOGSUM_SCALE) ];
}

/* Takes two log numbers and returns the log of their real sum (exact) */
// static inline
float calc_Logsum_exact(float x, float y)
{
   return log( exp(x) + exp(y) );
}

/* Print the logsum table */
void print_Logsum()
{
   printf("=== LOGSUM TABLE ===\n");
   for (int i = 0; i < 16000; i += 160)
   {
      printf("[%d] %.2f\n", i, logsum_lookup[i]);
   }
   printf("\n\n");
}

/* Get the number of characters in a string (including \0) */
int get_str_len(char* str) 
{
   int i;
   while (str[i] != '\0') {
      i++;
   } 
   i++;
   return i;
}

/*
 *  FUNCTION:  dp_matrix_Print()
 *  SYNOPSIS:  Print out dynamic programming matrix to screen.
 *
 *  PURPOSE:
 *
 *  ARGS:      <Q>         query length,
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix
 *
 *  RETURN:
 */
void dp_matrix_Print (const int Q, const int T,
                      const float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                      const float sp_MX[ NUM_SPECIAL_STATES * (Q + 1) ])
{
   /* PRINT resulting dp matrix */
   printf("\n\n==== DP MATRIX BEGIN ==== \n");
   /* Header */
   printf("\t");
   for (int i = 0; i <= T; i++)
   {
      printf("%d\t", i);
   }
   printf("\n");

   /* Row-by-Row */
   for (int i = 0; i < Q + 1; i++)
   {
      printf( "%d M\t", i );
      for (int j = 0; j <= T; j++)
      {
         printf( "%.3f\t", MMX(i, j) );
      }
      printf("\n");

      printf( "%d I\t", i );
      for (int j = 0; j <= T; j++)
      {
         printf( "%.3f\t", IMX(i, j) );
      }
      printf("\n");

      printf( "%d D\t", i );
      for (int j = 0; j <= T; j++)
      {
         printf( "%.3f\t", DMX(i, j) );
      }
      printf("\n\n");
   }

   printf("=== SPECIAL STATES ===\n");
   printf("N:\t");
   for (int i = 0; i <= Q; i++)
   { printf( "%.3f\t", XMX(SP_N, i) ); }
   printf("\n");
   printf("J:\t");
   for (int i = 0; i <= Q; i++)
   { printf( "%.3f\t", XMX(SP_J, i) ); }
   printf("\n");
   printf("E:\t");
   for (int i = 0; i <= Q; i++)
   { printf( "%.3f\t", XMX(SP_E, i) ); }
   printf("\n");
   printf("C:\t");
   for (int i = 0; i <= Q; i++)
   { printf( "%.3f\t", XMX(SP_C, i) ); }
   printf("\n");
   printf("B:\t");
   for (int i = 0; i <= Q; i++)
   { printf( "%.3f\t", XMX(SP_B, i) ); }
   printf("\n");

   printf("==== DP MATRIX END ==== \n\n");
}

/*
 *  FUNCTION:  test_matrix_Print()
 *  SYNOPSIS:  Print out dynamic programming matrix to screen.
 *
 *  PURPOSE:
 *
 *  ARGS:      <Q>         query length,
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix
 *
 *  RETURN:
 */
void test_matrix_Print (const int Q, const int T,
                        const float test_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ])
{
   /* PRINT resulting dp matrix */
   printf("\n\n==== TEST MATRIX BEGIN ==== \n");
   /* Header */
   printf("\t");
   for (int i = 0; i <= T; i++)
   {
      printf("%d\t", i);
   }
   printf("\n");

   /* Row-by-Row */
   for (int i = 0; i < Q + 1; i++)
   {
      printf( "%d\t", i );
      for (int j = 0; j <= T; j++)
      {
         printf( "%3.0f\t", TMX(i, j) );
      }
      printf("\n");
   }

   printf("==== TEST MATRIX END ==== \n\n");
}

/* Clear all matrix values to -INF. (for testing) */
void dp_matrix_Clear (const int Q, const int T,
                      float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                      float sp_MX[ NUM_SPECIAL_STATES * (Q + 1) ])
{
   for (int i = 0; i <= Q; i++)
   {
      for (int j = 0; j <= T; j++) {
         MMX(i, j) = IMX(i, j) = DMX(i, j) = -INF;
      }

      for (int j = 0; j < NUM_SPECIAL_STATES; j++) {
         XMX(j, i) = -INF;
      }
   }
}

/* Clear all matrix values to -INF. (for testing) */
void dp_matrix_Clear3 (const int Q, const int T,
                      float st_MX3[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                      float sp_MX[ NUM_SPECIAL_STATES * (Q + 1) ])
{
   for (int i = 0; i <= Q; i++)
   {
      for (int j = 0; j < NUM_SPECIAL_STATES; j++) {
         XMX(j, i) = -INF;
      }
   }

   for (int i = 0; i < 3; i++) {
      for (int j = 0; j <= T; j++) {
         MMX3(i, j) = IMX3(i, j) = DMX3(i, j) = -INF;
      }
   }
}

/* Set all matrix values to val */
void dp_matrix_Clear_X (const int Q, const int T,
                        float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                        float sp_MX[ NUM_SPECIAL_STATES * (Q + 1) ],
                        float val)
{
   for (int i = 0; i <= Q; i++)
   {
      for (int j = 0; j <= T; j++) {
         MMX(i, j) = IMX(i, j) = DMX(i, j) = val;
      }

      for (int j = 0; j < NUM_SPECIAL_STATES; j++) {
         XMX(j, i) = val;
      }
   }
}


/* Set all matrix values to val */
int dp_matrix_Compare (const int Q, const int T,
                        float st_MX_1[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                        float sp_MX_1[ NUM_SPECIAL_STATES * (Q + 1) ],
                        float st_MX_2[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                        float sp_MX_2[ NUM_SPECIAL_STATES * (Q + 1) ] )
{
   int i, j, st;

   for (i = 0; i <= Q; i++)
   {
      for (j = 0; j <= T; j++) 
      {
         for (st = 0; st < NUM_NORMAL_STATES; st++) 
         {
            if ( ST_MX(st_MX_1, st, i, j) != ST_MX(st_MX_1, st, i, j) ) 
            {
               return false;
            } 
         }
      }

      for (st = 0; st < NUM_SPECIAL_STATES; st++) 
      {
         if ( SP_MX(sp_MX_1, st, i) != SP_MX(sp_MX_2, st, i) ) 
         {
            return false;
         }
      }
   }
   return true;
}

/* Copy source matrix into destination matrix */
void dp_matrix_Copy (const int Q, const int T,
                     float st_MX_src[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                     float sp_MX_src[ NUM_SPECIAL_STATES * (Q + 1) ],
                     float st_MX_dst[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                     float sp_MX_dst[ NUM_SPECIAL_STATES * (Q + 1) ] )
{
   int i, j, st;

   for (i = 0; i <= Q; i++)
   {
      for (j = 0; j <= T; j++) 
      {
         for (st = 0; st < NUM_NORMAL_STATES; st++) 
         {
            // printf("(%d, %d, %d)\n", i, j, st);
            // printf("(%f)\n", st_MX_src[0]);
            // printf("(%f)\n", st_MX_dst[0]);

            ST_MX(st_MX_dst, st, i, j) = ST_MX(st_MX_src, st, i, j);
         }
      }

      for (st = 0; st < NUM_SPECIAL_STATES; st++) 
      {
         // printf("(%d, %d)\n", i, st);
         SP_MX(sp_MX_dst, st, i) = SP_MX(sp_MX_src, st, i);
      }
   }
}


/*
 *  FUNCTION:  dp_matrix_Save()
 *  SYNOPSIS:  Save dynamic programming matrix to file.
 *
 *  PURPOSE:
 *
 *  ARGS:      <Q>         query length,
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State
 *             <f>         Filename
 *
 *  RETURN:
 */
void dp_matrix_Save (const int Q, const int T, 
                           const float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                           const float sp_MX[ NUM_SPECIAL_STATES * (Q + 1) ],
                           const char *_filename_)
{
   FILE *fp;
   fp = fopen(_filename_, "w");

   const char* STATE_NAMES[] = {
   "M_ST",
   "I_ST",
   "D_ST",
   "E_ST",
   "N_ST",
   "J_ST",
   "C_ST",
   "B_ST",
   "S_ST",
   "T_ST",
   "X_ST",
   };

   /* PRINT resulting dp matrix */
   fprintf(fp, "##### DP MATRIX ##### \n");
   fprintf(fp, "XDIM\t%d\t%d\n\n", Q, T);

   /* Header */
   fprintf(fp, "##### NORMAL STATES #####\n");
   fprintf(fp, "XMATRIX\n");
   /* Header Indices */
   fprintf(fp, "#\t");
   for (int i = 0; i <= T; i++)
   {
      fprintf(fp, "%d\t", i);
   }
   fprintf(fp, "\n");

   /* Row-by-Row Values */
   for (int i = 0; i < Q + 1; i++)
   {
      fprintf(fp, "M %d\t", i );
      for (int j = 0; j <= T; j++)
      {
         fprintf(fp, "%.3f\t", MMX(i, j) );
      }
      fprintf(fp, "\n");

      fprintf(fp, "I %d\t", i );
      for (int j = 0; j <= T; j++)
      {
         fprintf(fp, "%.3f\t", IMX(i, j) );
      }
      fprintf(fp, "\n");

      fprintf(fp, "D %d\t", i );
      for (int j = 0; j <= T; j++)
      {
         fprintf(fp, "%.3f\t", DMX(i, j) );
      }
      fprintf(fp, "\n\n");
   }
   fprintf(fp, "/\n\n");

   fprintf(fp, "###### SPECIAL STATES #####\n");
   fprintf(fp, "N\t");
   for (int i = 0; i <= Q; i++)
   { fprintf(fp, "%.3f\t", XMX(SP_N, i) ); }
   fprintf(fp, "\n");
   fprintf(fp, "J\t");
   for (int i = 0; i <= Q; i++)
   { fprintf(fp, "%.3f\t", XMX(SP_J, i) ); }
   fprintf(fp, "\n");
   fprintf(fp, "E\t");
   for (int i = 0; i <= Q; i++)
   { fprintf(fp, "%.3f\t", XMX(SP_E, i) ); }
   fprintf(fp, "\n");
   fprintf(fp, "C\t");
   for (int i = 0; i <= Q; i++)
   { fprintf(fp, "%.3f\t", XMX(SP_C, i) ); }
   fprintf(fp, "\n");
   fprintf(fp, "B\t");
   for (int i = 0; i <= Q; i++)
   { fprintf(fp, "%.3f\t", XMX(SP_B, i) ); }
   fprintf(fp, "\n");

   fclose(fp);
}

/*
 *  FUNCTION:  dp_matrix_Save()
 *  SYNOPSIS:  Save dynamic programming matrix to file.
 *
 *  PURPOSE:
 *
 *  ARGS:      <Q>         query length,
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State
 *             <tr>        TRACE object
 *             <f>         Filename
 *
 *  RETURN:
 */
void dp_matrix_trace_Save (const int Q, const int T, 
                           const float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                           const float sp_MX[ NUM_SPECIAL_STATES * (Q + 1) ],
                           const TRACEBACK *tr,
                           const char *_filename_)
{
   printf("Saving matrix...\n");

   FILE *fp;
   fp = fopen(_filename_, "w");

   const char* STATE_NAMES[] = {
   "M_ST",
   "I_ST",
   "D_ST",
   "E_ST",
   "N_ST",
   "J_ST",
   "C_ST",
   "B_ST",
   "S_ST",
   "T_ST",
   "X_ST",
   };

   /* PRINT resulting dp matrix */
   fprintf(fp, "##### DP MATRIX ##### \n");
   fprintf(fp, "XDIM\t%d\t%d\n\n", Q, T);

   /* Traceback */
   fprintf(fp, "XTRACE\n");
   /* Traceback */
   for (int i = 0; i < tr->N; i++) 
   {
      int st = tr->traces[i].st;
      if ( st == M_ST ) {
         fprintf(fp, "%s\t%d\t%d\n", STATE_NAMES[st], tr->traces[i].i, tr->traces[i].j);
      }
   }
   fprintf(fp, "/\n\n");

   /* Header */
   fprintf(fp, "##### NORMAL STATES #####\n");
   fprintf(fp, "XMATRIX\n");
   /* Header Indices */
   fprintf(fp, "#\t");
   for (int i = 0; i <= T; i++)
   {
      fprintf(fp, "%d\t", i);
   }
   fprintf(fp, "\n");

   /* Row-by-Row Values */
   for (int i = 0; i <= Q; i++)
   {
      fprintf(fp, "M %d\t", i );
      for (int j = 0; j <= T; j++)
      {
         fprintf(fp, "%.3f\t", MMX(i, j) );
      }
      fprintf(fp, "\n");

      fprintf(fp, "I %d\t", i );
      for (int j = 0; j <= T; j++)
      {
         fprintf(fp, "%.3f\t", IMX(i, j) );
      }
      fprintf(fp, "\n");

      fprintf(fp, "D %d\t", i );
      for (int j = 0; j <= T; j++)
      {
         fprintf(fp, "%.3f\t", DMX(i, j) );
      }
      fprintf(fp, "\n\n");
   }
   fprintf(fp, "/\n\n");

   fprintf(fp, "###### SPECIAL STATES #####\n");
   fprintf(fp, "XSPECIAL\n");
   fprintf(fp, "N\t");
   for (int i = 0; i <= Q; i++)
   { 
      fprintf(fp, "%.3f\t", XMX(SP_N, i) ); 
   }
   fprintf(fp, "\n");
   fprintf(fp, "J\t");
   for (int i = 0; i <= Q; i++)
   { 
      fprintf(fp, "%.3f\t", XMX(SP_J, i) ); 
   }
   fprintf(fp, "\n");
   fprintf(fp, "E\t");
   for (int i = 0; i <= Q; i++)
   { 
      fprintf(fp, "%.3f\t", XMX(SP_E, i) ); 
   }
   fprintf(fp, "\n");
   fprintf(fp, "C\t");
   for (int i = 0; i <= Q; i++)
   {
      fprintf(fp, "%.3f\t", XMX(SP_C, i) ); 
   }
   fprintf(fp, "\n");
   fprintf(fp, "B\t");
   for (int i = 0; i <= Q; i++)
   { 
      fprintf(fp, "%.3f\t", XMX(SP_B, i) ); 
   }
   fprintf(fp, "\n/\n");

   fclose(fp);

   printf("Saved Matrix to: '%s'\n", _filename_);
}