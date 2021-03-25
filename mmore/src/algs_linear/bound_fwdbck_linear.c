/*******************************************************************************
 *  FILE:      bounded_fwdbck_linear.c
 *  PURPOSE:   Bounded Forward/Backward Algorithm 
 *             (Linear Space Alg)
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/* local imports */
#include "../objects/structs.h"
#include "../utilities/_utilities.h"
#include "../objects/_objects.h"

/* self header */
#include "_algs_linear.h"
#include "bound_fwdbck_linear.h"

/*
 *      NOTE: HOW TO CONVERT row-coords to diag-coords
 *       MMX3(i-1,j-1) => MMX3(, d_2)
 *       MMX3(i,  j-1) => MMX3(, d_1)
 *       MMX3(i,  j  ) => MMX3(, d_1)
 */

/* private functions */
static inline 
float 
MY_Sum( const float x, const float y );

static inline 
float 
MY_Prod( const float x, const float y );

static inline 
float 
MY_Zero();

static inline 
float 
MY_One();

/*! FUNCTION: run_Bound_Forward_Linear()
 *  SYNOPSIS: Perform Edge-Bounded Forward step of Cloud Search Algorithm.
 *            Runs traditional Forward-Backward Algorithm, but only performs
 *             computation on cells that fall within the bounds determined by
 *             the <edg> EDGEBOUNDS object, which stores a series of 
 *             (left-bound, right-bound) pairs sorted by row.  
 *            Normal state matrix is stored in linear space.
 *             <st_MX3> is size [3 * (Q + T + 1)]. Only requires size [2 * (T + 1)],
 *             but is reused from cloud_forward_().
 *            Final score produced by Forward is stored in <sc_final>.
 *
 *  RETURN:   Returns the final score of the Forward Algorithm.
 */
int 
run_Bound_Forward_Linear(     const SEQUENCE*            query,         /* query sequence */
                              const HMM_PROFILE*         target,        /* target HMM model */
                              const int                  Q,             /* query length */
                              const int                  T,             /* target length */
                              MATRIX_3D* restrict        st_MX3,        /* normal state matrix */
                              MATRIX_2D* restrict        sp_MX,         /* special state matrix */
                              EDGEBOUNDS*                edg,           /* edgebounds */
                              float*                     sc_final )     /* (OUTPUT) final score */
{
   // /* Local Method for computing sum */
   // float (*MY_Sum)(const float x, const float y) = MATH_LogSum;
   // /* Local Method for computing product */
   // float (*MY_Prod)(const float x, const float y) = MATH_LogProd;
   // /* Local Method for computing one */
   // float (*MY_One)() = MATH_LogOne;
   // /* Local Method for computing zero */
   // float (*MY_Zero)() = MATH_LogZero;

   /* vars for accessing query/target data structs */
   char     a;                               /* store current character in sequence */
   int      A;                               /* store int value of character */
   char*    seq;                             /* alias for getting seq */
   int      N;                               /* length of edgebound list */
   bool     is_local;                        /* whether using local or global alignments */

   /* vars for indexing into data matrices by row-col */
   int      b, d, i, j, k;                   /* antidiagonal, row, column indices */
   int      q_0, q_1;                        /* real index of current and previous rows (query) */
   int      qx0, qx1;                        /* maps column index into data index (query) */
   int      t_0, t_1;                        /* real index of current and previous columns (target) */
   int      tx0, tx1;                        /* maps target index into data index (target)  */
   int      t_range;                         /* range of targets on current row */

   /* vars for indexing into data matrices by anti-diag */
   int      d_0, d_1, d_2;                   /* real index of current and previous antidiagonals */
   int      dx0, dx1, dx2;                   /* mod mapping of antidiagonal index into data matrix */
   int      k_0, k_1;                        /* offset into antidiagonal */
   int      d_st, d_end, d_cnt;              /* starting and ending diagonal indices */
   int      dim_T, dim_Q, dim_TOT;           /* dimensions of submatrix being searched */
   int      dim_min, dim_max;                /* diagonal index where num cells reaches highest point and diminishing point */ 
   int      num_cells;                       /* number of cells in current diagonal */

   /* vars for indexing into edgebound lists */
   BOUND    bnd;                             /* current bound */
   int      id;                              /* id in edgebound list (row/diag) */
   int      r_0;                             /* current index for current row */
   int      r_0b, r_0e;                      /* begin and end indices for current row in edgebound list */
   int      r_1;                             /* current index for previous row */
   int      r_1b, r_1e;                      /* begin and end indices for current row in edgebound list */
   int      le_0, re_0;                      /* right/left matrix bounds of current diag */
   int      lb_0, rb_0;                      /* bounds of current search space on current diag */
   int      lb_1, rb_1;                      /* bounds of current search space on previous diag */
   int      lb_2, rb_2;                      /* bounds of current search space on 2-back diag */
   bool     rb_T;                            /* checks if edge touches right bound of matrix */

   /* vars for recurrance scores */
   float    prv_M, prv_I, prv_D;             /* previous (M) match, (I) insert, (D) delete states */
   float    prv_B, prv_E;                    /* previous (B) begin and (E) end states */
   float    prv_J, prv_N, prv_C;             /* previous (J) jump, (N) initial, and (C) terminal states */
   float    prv_loop, prv_move;              /* previous loop and move for special states */
   float    prv_sum, prv_best;               /* temp subtotaling vars */
   float    sc_best;                         /* final best scores */
   float    sc_M, sc_I, sc_D, sc_E;          /* match, insert, delete, end scores */

   /* for domain ranges */
   RANGE*               dom_range;           /* domain range */          
   RANGE                T_range;             /* target range */
   RANGE                Q_range;             /* query range */
   bool                 is_q_0_in_dom_range; /* checks if current query position is inside the domain range */
   bool                 is_q_1_in_dom_range; /* checks if previous query position is inside the domain range */

   /* debugger tools */
   FILE*       dbfp;
   MATRIX_2D*  cloud_MX;
   MATRIX_2D*  cloud_MX3;
   MATRIX_3D*  test_MX;
   MATRIX_3D*  test_MX3;
   EDGEBOUNDS* test_edg;
   int         num_writes;
   int         num_clears;

   /* initialize debugging matrix */
   #if DEBUG
   {
      cloud_MX    = debugger->cloud_MX;
      cloud_MX3   = debugger->cloud_MX3;
      test_MX     = debugger->test_MX;
      test_MX3    = debugger->test_MX3;
      test_edg    = debugger->test_edg;

      MATRIX_2D_Reuse( cloud_MX, Q+1, T+1 );
      MATRIX_2D_Fill( cloud_MX, 0 );
      MATRIX_2D_Reuse( cloud_MX3, 3, (Q+1)+(T+1) );
      MATRIX_2D_Fill( cloud_MX3, 0 );
      MATRIX_3D_Reuse( test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
      MATRIX_3D_Fill( test_MX, -INF );
      MATRIX_3D_Reuse( test_MX3, NUM_NORMAL_STATES, 3, (Q+1)+(T+1) );
      MATRIX_3D_Fill( test_MX3, -INF );
      EDGEBOUNDS_Reuse( test_edg, Q, T );

      num_writes = 0;
      num_clears = 0;
   }
   #endif

   #if DEBUG 
   {
      int dirty_vals = MATRIX_3D_Check_Clean( st_MX3 );
      printf("PRE-CLEAN CHECK -> BOUND FWD LIN:\t%d\n", dirty_vals);
      if (dirty_vals > 0) {
         MATRIX_3D_Dump( st_MX3, stdout );
      }
   }
   #endif

   /* --------------------------------------------------------------------------------- */

   /* initialize logsum lookup table if it has not already been */
   MATH_Logsum_Init();

   /* query sequence */
   seq         = query->seq;
   N           = EDGEBOUNDS_GetSize( edg );
   /* local or global alignments? */
   is_local    = target->isLocal;
   sc_E        = (is_local) ? 0 : -INF;
   if ( is_local == true ) {
      sc_E = MY_One();
   } else {
      sc_E = MY_Zero();
   }

   /* domain range (query sequence) */
   Q_range.beg = 0;
   Q_range.end = Q + 1;
   /* target range */
   T_range.beg = 0;
   T_range.end = T + 1;

   /* UNROLLED INTIAL ROW */
   q_0 = Q_range.beg;
   {
      qx0 = q_0 % 2;

      /* get edgebound range */
      r_0b  = EDGEBOUNDS_GetIndex_byRow_Fwd( edg, q_0 );
      r_0e  = EDGEBOUNDS_GetIndex_byRow_Fwd( edg, q_0 + 1 );

      /* initialize 0 row (top-edge) */
      for (t_0 = 0; t_0 < T; t_0++) 
      { 
         MMX3(qx0, t_0) = MY_Zero();
         IMX3(qx0, t_0) = MY_Zero();
         DMX3(qx0, t_0) = MY_Zero();
      }

      /* Initialize special states (?) */
      XMX(SP_N, q_0) = MY_One();                /* S->N, p=1             */
      XMX(SP_B, q_0) = XSC(SP_N,SP_MOVE);       /* S->N->B, no N-tail    */
      XMX(SP_E, q_0) = MY_Zero();
      XMX(SP_C, q_0) = MY_Zero();
      XMX(SP_J, q_0) = MY_Zero();               /* need seq to get here (?)  */

      /* init lookback 1 row */
      r_1b = r_0b;
      r_1e = r_0e;
   }

   /* MAIN RECURSION */
   /* FOR every position in QUERY sequence (row in matrix) */
   for (q_0 = Q_range.beg + 1; q_0 < Q_range.end; q_0++)
   {
      q_1 = q_0 - 1;
      qx0 = q_0 % 2;
      qx1 = q_1 % 2;

      t_0 = 0;
      tx0 = t_0;

      /* get edgebound range */
      r_0b  = EDGEBOUNDS_GetIndex_byRow_Fwd( edg, q_0 );
      r_0e  = EDGEBOUNDS_GetIndex_byRow_Fwd( edg, q_0 + 1 );

      /* Get next sequence character */
      a = seq[q_1];  /* off-by-one */
      A = AA_REV[a];

      /* Initialize E state for row */
      XMX(SP_E, q_0) = MY_Zero();

      /* FOR every EDGEBOUND in current ROW */
      for (r_0 = r_0b; r_0 < r_0e; r_0++)
      {
         /* in this context, "id" represents the "row" */
         bnd   = EDGEBOUNDS_Get(edg, r_0);
         id    = bnd.id;                           /* NOTE: this is always the same as cur_row, q_0 */
         lb_0  = MAX(bnd.lb - 1, T_range.beg);     /* can't overflow the left edge */
         rb_0  = MIN(bnd.rb, T_range.end);         /* can't overflow the right edge */

         /* FOR first position in EDGEBOUND profile */
         t_0 = lb_0;
         {
            t_1 = t_0 - 1;
            tx0 = t_0;
            tx1 = tx0 - 1;

            /* Initialize zero column (left-edge) */
            MMX3(qx0, tx0) = MY_Zero();
            IMX3(qx0, tx0) = MY_Zero();
            DMX3(qx0, tx0) = MY_Zero();

            /* embed linear row into quadratic test matrix */
            #if DEBUG
            {
               MX_2D(cloud_MX, q_0, t_0) = 1.0;
               MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX3(qx0, t_0);
               MX_3D(test_MX, INS_ST, q_0, t_0) = IMX3(qx0, t_0);
               MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX3(qx0, t_0);
            }
            #endif
         }

         /* MAIN RECURSION */
         /* FOR every position in TARGET profile */
         for (t_0 = lb_0 + 1; t_0 < rb_0 - 1; t_0++)
         {
            t_1 = t_0 - 1;
            tx0 = t_0;
            tx1 = tx0 - 1;

            /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
            /* best previous state transition (match takes the diag element of each prev state) */
            prv_M = MY_Prod( MMX3(qx1, tx1), TSC(t_1, M2M) );
            prv_I = MY_Prod( IMX3(qx1, tx1), TSC(t_1, I2M) );
            prv_D = MY_Prod( DMX3(qx1, tx1), TSC(t_1, D2M) );
            prv_B = MY_Prod( XMX(SP_B, q_1), TSC(t_1, B2M) ); /* from begin match state (new alignment) */
            /* best-to-match */
            prv_sum = MY_Sum( MY_Sum( prv_M, prv_I ),
                              MY_Sum( prv_B, prv_D ) );
            MMX3(qx0, tx0) = MY_Prod( prv_sum, MSC(t_0, A) );

            /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
            /* previous states (match takes the previous row (upper) of each state) */
            prv_M = MY_Prod( MMX3(qx1, tx0), TSC(t_0, M2I) );
            prv_I = MY_Prod( IMX3(qx1, tx0), TSC(t_0, I2I) );
            /* best-to-insert */
            prv_sum = MY_Sum( prv_M, prv_I );
            IMX3(qx0, tx0) = MY_Prod( prv_sum, ISC(t_0, A) );

            /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) */
            /* previous states (match takes the previous column (left) of each state) */
            prv_M = MY_Prod( MMX3(qx0, tx1), TSC(t_1, M2D) );
            prv_D = MY_Prod( DMX3(qx0, tx1), TSC(t_1, D2D) );
            /* best-to-delete */
            prv_sum = MY_Sum( prv_M, prv_D );
            DMX3(qx0, tx0) = prv_sum;

            /* UPDATE E STATE */
            prv_M = MY_Prod( MMX3(qx0, tx0), sc_E );
            prv_D = MY_Prod( DMX3(qx0, tx0), sc_E );
            /* best-to-e-state */
            prv_E = XMX(SP_E, q_0);
            XMX(SP_E, q_0) = MY_Sum( MY_Sum( prv_M, prv_D ),
                                             prv_E );

            /* embed linear row into quadratic test matrix */
            #if DEBUG
            {
               MX_2D(cloud_MX, q_0, t_0) = 1.0;
               MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX3(qx0, t_0);
               MX_3D(test_MX, INS_ST, q_0, t_0) = IMX3(qx0, t_0);
               MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX3(qx0, t_0);
            }
            #endif
         }

         /* UNROLLED FINAL LOOP ITERATION */
         if ( rb_0 - lb_0 > 1 )
         {
            t_0   = rb_0 - 1;
            t_1   = t_0 - 1;
            tx0   = t_0;
            tx1   = tx0 - 1;

            /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
            /* best previous state transition (match takes the diag element of each prev state) */
            prv_M = MY_Prod( MMX3(qx1, tx1), TSC(t_1, M2M) );
            prv_I = MY_Prod( IMX3(qx1, tx1), TSC(t_1, I2M) );
            prv_D = MY_Prod( DMX3(qx1, tx1), TSC(t_1, D2M) );
            prv_B = MY_Prod( XMX(SP_B, q_1), TSC(t_1, B2M) );    /* from begin match state (new alignment) */
            /* sum-to-match */
            prv_sum = MY_Sum( MY_Sum( prv_M, prv_I ),
                              MY_Sum( prv_D, prv_B ) );
            MMX3(qx0, tx0) = prv_sum + MSC(t_0, A);

            /* FIND SUM OF PATHS TO INSERT STATE (unrolled) */
            IMX3(qx0, tx0) = MY_Zero();

            /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) (unrolled) */
            /* previous states (match takes the left element of each state) */
            prv_M = MY_Prod( MMX3(qx0, tx1), TSC(t_1, M2D) );
            prv_D = MY_Prod( DMX3(qx0, tx1), TSC(t_1, D2D) );
            /* sum-to-delete */
            prv_sum = MY_Sum( prv_M, prv_D );
            DMX3(qx0, tx0) = prv_sum;

            /* UPDATE E STATE (unrolled) */
            prv_E = XMX(SP_E, q_0);
            prv_M = MMX3(qx0, t_0);
            prv_D = DMX3(qx0, t_0);
            /* best-to-begin */
            XMX(SP_E, q_0) = MY_Sum( MY_Sum( prv_D, prv_M ),
                                             prv_E ); 

            /* embed linear row into quadratic test matrix */
            #if DEBUG
            {
               MX_2D(cloud_MX, q_0, t_0) = 1.0;
               MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX3(qx0, t_0);
               MX_3D(test_MX, INS_ST, q_0, t_0) = IMX3(qx0, t_0);
               MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX3(qx0, t_0);
            }
            #endif
         }
      }

      /* SPECIAL STATES */
      /* J state */
      prv_J = MY_Prod( XMX(SP_J, q_1), XSC(SP_J, SP_LOOP) );       /* J->J */
      prv_E = MY_Prod( XMX(SP_E, q_0), XSC(SP_E, SP_LOOP) );       /* E->J is E's "loop" */
      XMX(SP_J, q_0) = MY_Sum( prv_J, prv_E );         

      /* C state */
      prv_C = MY_Prod( XMX(SP_C, q_1), XSC(SP_C, SP_LOOP) );
      prv_E = MY_Prod( XMX(SP_E, q_0), XSC(SP_E, SP_MOVE) );
      XMX(SP_C, q_0) = MY_Sum( prv_C, prv_E );

      /* N state */
      prv_N = MY_Prod( XMX(SP_N, q_1), XSC(SP_N, SP_LOOP) );
      XMX(SP_N, q_0) = prv_N;

      /* B state */
      prv_N = MY_Prod( XMX(SP_N, q_0), XSC(SP_N, SP_MOVE) );         /* N->B is N's move */
      prv_J = MY_Prod( XMX(SP_J, q_0), XSC(SP_J, SP_MOVE) );         /* J->B is J's move */
      XMX(SP_B, q_0) = MY_Sum( prv_N, prv_J );    

      /* SCRUB PREVIOUS ROW */
      for (r_1 = r_1b; r_1 < r_1e; r_1++) 
      {
         /* in this context, "id" represents the "row" */
         bnd   = EDGEBOUNDS_Get(edg, r_1);      /* NOTE: this is always the same as cur_row, q_0 */
         lb_1  = MAX(bnd.lb - 1, T_range.beg);     /* can't overflow the left edge */
         rb_1  = MIN(bnd.rb, T_range.end);         /* can't overflow the right edge */

         for (t_0 = lb_1; t_0 < rb_1; t_0++) {
            MMX3(qx1, t_0) = MY_Zero();
            IMX3(qx1, t_0) = MY_Zero();
            DMX3(qx1, t_0) = MY_Zero();
         }
      }

      /* SET CURRENT ROW TO PREVIOUS ROW */
      r_1b = r_0b;
      r_1e = r_0e;
   }

   /* UNROLLED FINAL LOOP (just cleanup) */
   q_0 = Q;
   {
      q_1 = Q - 1;
      qx0 = Q % 2;
      qx1 = Q % 2;

      /* SCRUB FINAL ROW */
      for (r_1 = r_1b; r_1 < r_1e; r_1++) 
      {
         /* in this context, "id" represents the "row" */
         bnd   = EDGEBOUNDS_Get(edg, r_1);        /* NOTE: this is always the same as cur_row, q_0 */
         lb_1  = MAX(bnd.lb - 1, T_range.beg);     /* can't overflow the left edge */
         rb_1  = MIN(bnd.rb, T_range.end);         /* can't overflow the right edge */

         for (t_0 = lb_1; t_0 < rb_1; t_0++) {
            MMX3(qx1, t_0) = MY_Zero();
            IMX3(qx1, t_0) = MY_Zero();
            DMX3(qx1, t_0) = MY_Zero();
         }
      }
   }
   
   /* T state */
   sc_best     = MY_Prod( XMX(SP_C, Q), XSC(SP_C, SP_MOVE) );
   *sc_final   = sc_best;

   /* flag matrices that they contain dirty values (not -INF) */
   st_MX3->clean  = true;
   sp_MX->clean   = true;

   #if DEBUG 
   {
      int dirty_vals = MATRIX_3D_Check_Clean( st_MX3 );
      printf("POST-CLEAN CHECK -> BOUND FWD LIN:\t%d\n", dirty_vals);
      // if (dirty_vals > 0) {
      //    MATRIX_3D_Dump( st_MX3, stdout );
      // }
   }
   #endif

   return STATUS_SUCCESS;
}


/*! FUNCTION: run_Bound_Backward_Linear()
 *  SYNOPSIS: Perform Edge-Bounded Backward step of Cloud Search Algorithm.
 *            Runs traditional Forward-Backward Algorithm, but only performs
 *             computation on cells that fall within the bounds determined by
 *             the <edg> EDGEBOUNDS object, which stores a series of 
 *             (left-bound, right-bound) pairs sorted by row.  
 *            Normal state matrix is stored in linear space.
 *             <st_MX3> is size [3 * (Q + T + 1)]. Only requires size [2 * (T + 1)],
 *             but is reused from cloud_forward_().
 *            Final score produced by Backward is stored in <sc_final>.
 *
 *  RETURN:   Returns the final score of the Backward Algorithm.
 */
int run_Bound_Backward_Linear(   const SEQUENCE*         query,         /* query sequence */
                                 const HMM_PROFILE*      target,        /* target HMM model */
                                 const int               Q,             /* query length */
                                 const int               T,             /* target length */
                                 MATRIX_3D* restrict     st_MX3,        /* normal state matrix */
                                 MATRIX_2D* restrict     sp_MX,         /* special state matrix */
                                 EDGEBOUNDS*             edg,           /* edgebounds */
                                 float*                  sc_final )     /* (OUTPUT) final score */
{
   // /* Local Method for computing sum */
   // float (*MY_Sum)(const float x, const float y) = MATH_LogSum;
   // /* Local Method for computing product */
   // float (*MY_Prod)(const float x, const float y) = MATH_LogProd;
   // /* Local Method for computing one */
   // float (*MY_One)() = MATH_LogOne;
   // /* Local Method for computing zero */
   // float (*MY_Zero)() = MATH_LogZero;

   /* vars for accessing query/target data structs */
   char     a;                               /* store current character in sequence */
   int      A;                               /* store int value of character */
   char*    seq;                             /* alias for getting seq */
   int      N;                               /* length of edgebound list */
   bool     is_local;                        /* whether using local or global alignments */

   /* vars for indexing into data matrices by row-col */
   int      b, d, i, j, k;                   /* antidiagonal, row, column indices */
   int      q_0, q_1;                        /* real index of current and previous rows (query) */
   int      qx0, qx1;                        /* mod mapping of column index into data matrix (query) */
   int      t_0, t_1;                        /* real index of current and previous columns (target) */
   int      tx0, tx1;                        /* maps target index into data index (target)  */

   /* vars for indexing into data matrices by anti-diag */
   int      d_0, d_1, d_2;                   /* real index of current and previous antidiagonals */
   int      dx0, dx1, dx2;                   /* mod mapping of antidiagonal index into data matrix */
   int      k_0, k_1;                        /* offset into antidiagonal */
   int      d_st, d_end, d_cnt;              /* starting and ending diagonal indices */
   int      dim_T, dim_Q, dim_TOT;           /* dimensions of submatrix being searched */
   int      dim_min, dim_max;                /* diagonal index where num cells reaches highest point and diminishing point */ 
   int      num_cells;                       /* number of cells in current diagonal */

   /* vars for indexing into edgebound lists */
   BOUND    bnd;                             /* current bound */
   int      id;                              /* id in edgebound list (row/diag) */
   int      r_0;                             /* current index for current row */
   int      r_0b, r_0e;                      /* begin and end indices for current row in edgebound list */
   int      r_1;                             /* current index for previous row */
   int      r_1b, r_1e;                      /* begin and end indices for current row in edgebound list */
   int      le_0, re_0;                      /* right/left matrix bounds of current diag */
   int      lb_0, rb_0;                      /* bounds of current search space on current diag */
   int      lb_1, rb_1;                      /* bounds of current search space on previous diag */
   int      lb_2, rb_2;                      /* bounds of current search space on 2-back diag */
   int      lb_T, rb_T;                      /* truncated bounds within domain */                           /* checks if edge touches right bound of matrix */

   /* vars for recurrance scores */
   float    prv_M, prv_I, prv_D;             /* previous (M) match, (I) insert, (D) delete states */
   float    prv_B, prv_E;                    /* previous (B) begin and (E) end states */
   float    prv_J, prv_N, prv_C;             /* previous (J) jump, (N) initial, and (C) terminal states */
   float    prv_loop, prv_move;              /* previous loop and move for special states */
   float    prv_sum, prv_best;               /* temp subtotaling vars */
   float    sc_best;                         /* final best scores */
   float    sc_M, sc_I, sc_D, sc_E;          /* match, insert, delete, end scores */

   /* for domain ranges */
   RANGE*               dom_range;           /* domain range */          
   RANGE                T_range;             /* target range */
   RANGE                Q_range;             /* query range */
   bool                 is_q_0_in_dom_range; /* checks if current query position is inside the domain range */
   bool                 is_q_1_in_dom_range; /* checks if previous query position is inside the domain range */

   /* debugger tools */
   FILE*       dbfp;
   MATRIX_2D*  cloud_MX;
   MATRIX_2D*  cloud_MX3;
   MATRIX_3D*  test_MX;
   MATRIX_3D*  test_MX3;
   EDGEBOUNDS* test_edg;
   int         num_writes;
   int         num_clears;

   /* initialize debugging matrix */
   #if DEBUG
   {
      cloud_MX    = debugger->cloud_MX;
      cloud_MX3   = debugger->cloud_MX3;
      test_MX     = debugger->test_MX;
      test_MX3    = debugger->test_MX3;
      test_edg    = debugger->test_edg;

      MATRIX_2D_Reuse( cloud_MX, Q+1, T+1 );
      MATRIX_2D_Fill( cloud_MX, 0 );
      MATRIX_2D_Reuse( cloud_MX3, 3, (Q+1)+(T+1) );
      MATRIX_2D_Fill( cloud_MX3, 0 );
      MATRIX_3D_Reuse( test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
      MATRIX_3D_Fill( test_MX, -INF );
      MATRIX_3D_Reuse( test_MX3, NUM_NORMAL_STATES, 3, (Q+1)+(T+1) );
      MATRIX_3D_Fill( test_MX3, -INF );
      EDGEBOUNDS_Reuse( test_edg, Q, T );

      num_writes = 0;
      num_clears = 0;
   }
   #endif

   #if DEBUG 
   {
      int dirty_vals = MATRIX_3D_Check_Clean( st_MX3 );
      printf("PRE-CLEAN CHECK -> BOUND BCK LIN:\t%d\n", dirty_vals);
      if (dirty_vals > 0) {
         MATRIX_3D_Dump( st_MX3, stdout );
      }
   }
   #endif

   /* --------------------------------------------------------------------------------- */

   /* initialize logsum lookup table if it has not already been */
   MATH_Logsum_Init();

   /* query sequence */
   seq         = query->seq;
   N           = EDGEBOUNDS_GetSize( edg );
   /* local or global alignments? */
   is_local    = target->isLocal;
   sc_E        = (is_local) ? 0 : -INF;
   if ( is_local == true ) {
      sc_E = MY_One();
   } else {
      sc_E = MY_Zero();
   }

   /* domain range (query sequence) */
   Q_range.beg = 0;
   Q_range.end = Q;
   /* target range */
   T_range.beg = 1;
   T_range.end = T;

   /* UNROLLED INITIAL ROW */
   q_0 = Q_range.end;
   {
      qx0 = q_0 % 2;

      /* get edgebound range */
      r_0b  = EDGEBOUNDS_GetIndex_byRow_Bck( edg, q_0 + 1 );
      r_0e  = EDGEBOUNDS_GetIndex_byRow_Bck( edg, q_0 );

      /* Initialize special states */
      XMX(SP_J, q_0) = MY_Zero();
      XMX(SP_B, q_0) = MY_Zero();
      XMX(SP_N, q_0) = MY_Zero();
      XMX(SP_C, q_0) = XSC(SP_C, SP_MOVE);
      XMX(SP_E, q_0) = MY_Prod( XMX(SP_C, q_0), XSC(SP_E, SP_MOVE) );

      /* if Q-row bounds are not empty and the right-most bound spans T ( covers bottom-right corner ) */
      bnd   = EDGEBOUNDS_Get( edg, r_0b );
      if ( (r_0b - r_0e > 0) && (bnd.rb > T_range.end) )
      {
         t_0   = T_range.end;
         tx0   = t_0;

         MMX3(qx0, tx0) = XMX(SP_E, q_0);
         IMX3(qx0, tx0) = MY_Zero();
         DMX3(qx0, tx0) = XMX(SP_E, q_0);

         #if DEBUG 
         {
            MX_2D(cloud_MX, q_0, t_0) += 2.0;
            MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX3(qx0, tx0);
            MX_3D(test_MX, INS_ST, q_0, t_0) = IMX3(qx0, tx0);
            MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX3(qx0, tx0);
         }
         #endif
      }

      /* Initialize normal states (sparse) */
      for (r_0 = r_0b; r_0 > r_0e; r_0--) 
      {
         /* get bound data */
         bnd =  EDGEBOUNDS_Get(edg, r_0);          /* bounds for current bound */
         lb_0 = MAX(bnd.lb, T_range.beg);          /* can't overflow the left edge */
         rb_0 = MIN(bnd.rb, T_range.end);          /* can't overflow the right edge */

         for (t_0 = rb_0-1; t_0 >= lb_0; t_0--)
         {
            /* real target index */
            t_1 = t_0 + 1;
            tx0 = t_0;
            tx1 = tx0 - 1;

            prv_E = MY_Prod( XMX(SP_E, Q), sc_E );
            prv_D = MY_Prod( DMX3(qx0, tx1), TSC(t_0, M2D) );
            MMX3(qx0, tx0) = MY_Sum( prv_E, prv_D );

            prv_E = MY_Prod( XMX(SP_E, Q), sc_E );
            prv_D = MY_Prod( DMX3(qx0, tx1), TSC(t_0, D2D) );
            DMX3(qx0, tx0) = MY_Sum( prv_E, prv_D );

            IMX3(qx0, tx0) = MY_Zero();

            #if DEBUG 
            {
               MX_2D(cloud_MX, q_0, t_0) += 2.0;
               MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX3(qx0, tx0);
               MX_3D(test_MX, INS_ST, q_0, t_0) = IMX3(qx0, tx0);
               MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX3(qx0, tx0);
            }  
            #endif
         }
      }

      /* init lookback 1 row */
      r_1b = r_0b;
      r_1e = r_0e;
   }

   /* MAIN RECURSION */
   /* FOR every bound in EDGEBOUND */
   for (q_0 = Q_range.end - 1; q_0 > Q_range.beg; q_0--)
   {
      q_1 = q_0 + 1;
      qx0 = q_0 % 2;
      qx1 = qx0 + 1;

      /* get edgebound range */
      r_0b  = EDGEBOUNDS_GetIndex_byRow_Bck( edg, q_0 + 1 );
      r_0e  = EDGEBOUNDS_GetIndex_byRow_Bck( edg, q_0 );

      /* Get next sequence character */
      a = seq[q_0];
      A = AA_REV[a];

      /* UPDATE SPECIAL STATES at the start of EACH ROW */

      /* B STATE (sparse) */
      XMX(SP_B, q_0) = MY_Zero();
      for (r_1 = r_1b; r_1 > r_1e; r_1--) 
      {
         bnd   = EDGEBOUNDS_Get(edg, r_1);          /* bounds for current bound */
         lb_0  = MAX(bnd.lb, T_range.beg);          /* can't overflow the left edge */
         rb_0  = MIN(bnd.rb, T_range.end);          /* can't overflow the right edge */

         for (t_0 = lb_0; t_0 <= rb_0; t_0++)
         {
            t_1   = t_0 - 1;
            tx0   = t_0;
            tx1   = tx0 - 1;

            prv_sum  = XMX(SP_B, q_0);
            prv_M    = MY_Prod( MY_Prod(  MMX3(qx1, tx0), TSC(t_1, B2M) ), 
                                          MSC(t_0, A) );
            XMX(SP_B, q_0) = MY_Sum( prv_sum, prv_M);
         }
      }

      prv_J = MY_Prod( XMX(SP_J, q_1), XSC(SP_J, SP_LOOP) );
      prv_B = MY_Prod( XMX(SP_B, q_0), XSC(SP_J, SP_MOVE) );
      XMX(SP_J, q_0) = MY_Sum( prv_J, prv_B );

      prv_C = MY_Prod( XMX(SP_C, q_1), XSC(SP_C, SP_LOOP) );
      XMX(SP_C, q_0) = prv_C;

      prv_J = MY_Prod( XMX(SP_J, q_0), XSC(SP_E, SP_LOOP) );
      prv_C = MY_Prod( XMX(SP_C, q_0), XSC(SP_E, SP_MOVE) );
      XMX(SP_E, q_0) = MY_Sum( prv_J, prv_C );

      prv_N = MY_Prod( XMX(SP_N, q_1), XSC(SP_N, SP_LOOP) );
      prv_B = MY_Prod( XMX(SP_B, q_0), XSC(SP_N, SP_MOVE) );
      XMX(SP_N, q_0) = MY_Sum( prv_N, prv_B );

      /* if there is a bound on row and the right-most bound spans T */
      bnd   = EDGEBOUNDS_Get( edg, r_0b );
      if ( (r_0b - r_0e > 0) && (bnd.rb > T_range.end) )
      {
         t_0 = T_range.end;
         tx0 = t_0;

         MMX3(qx0, tx0) = XMX(SP_E, q_0);
         IMX3(qx0, tx0) = MY_Zero();
         DMX3(qx0, tx0) = XMX(SP_E, q_0);

         #if DEBUG 
         {
            MX_2D(cloud_MX, q_0, t_0) += 2.0;
            MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX3(qx0, tx0);
            MX_3D(test_MX, INS_ST, q_0, t_0) = IMX3(qx0, tx0);
            MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX3(qx0, tx0);
         }
         #endif
      }

      /* FOR every EDGEBOUND in current ROW */
      for (r_0 = r_0b; r_0 > r_0e; r_0--)
      {
         /* get bound data */
         bnd   = EDGEBOUNDS_Get(edg, r_0);          /* bounds for current bound */
         lb_0  = MAX(bnd.lb, T_range.beg);          /* can't overflow the left edge */
         rb_0  = MIN(bnd.rb, T_range.end);          /* can't overflow the right edge */

         /* FOR every position in TARGET profile */
         for (t_0 = rb_0-1; t_0 >= lb_0; t_0--)
         {
            t_1   = t_0 + 1;
            tx0   = t_0;
            tx1   = tx0 - 1;

            /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
            prv_M = MY_Prod( MY_Prod(  MMX3(qx1, tx1), TSC(t_0, M2M) ),
                                       MSC(t_1, A) );
            prv_I = MY_Prod( MY_Prod(  IMX3(qx1, tx0), TSC(t_0, M2I) ),
                                       ISC(t_0, A) );
            prv_D = MY_Prod( DMX3(qx0, tx1), TSC(t_0, M2D) );
            prv_E = MY_Prod( XMX(SP_E, q_0), sc_E );     /* from end match state (new alignment) */
            /* best-to-match */
            prv_sum = MY_Sum(  MY_Sum( prv_M, prv_I ),
                                 MY_Sum( prv_E, prv_D ) );
            MMX3(qx0, tx0) = prv_sum;

            /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
            prv_M = MY_Prod( MY_Prod( MMX3(qx1, tx1), TSC(t_0, I2M) ),
                                      MSC(t_1, A) );
            prv_I = MY_Prod( MY_Prod( IMX3(qx1, tx0), TSC(t_0, I2I) ),
                                      ISC(t_0, A) );
            /* best-to-insert */
            prv_sum = MY_Sum( prv_M, prv_I );
            IMX3(qx0, tx0) = prv_sum;

            /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
            prv_M = MY_Prod( MY_Prod( MMX3(qx1, t_1), TSC(t_0, D2M)  ),
                                      MSC(t_1, A) );
            prv_D = DMX3(qx0, t_1) + TSC(t_0, D2D);
            prv_E = XMX(SP_E, q_0) + sc_E;
            /* best-to-delete */
            prv_sum = MY_Sum( prv_M, 
                      MY_Sum( prv_D, prv_E ) );
            DMX3(qx0, t_0) = prv_sum;

            #if DEBUG 
            {
               MX_2D(cloud_MX, q_0, t_0) += 2.0;
               MX_3D(test_MX, MAT_ST, q_0, t_0) = MMX3(qx0, t_0);
               MX_3D(test_MX, INS_ST, q_0, t_0) = IMX3(qx0, t_0);
               MX_3D(test_MX, DEL_ST, q_0, t_0) = DMX3(qx0, t_0);
            }
            #endif
         }
      }

      /* SCRUB PREVIOUS ROW */
      for (r_1 = r_1b; r_1 > r_1e; r_1--) 
      {
         bnd   = EDGEBOUNDS_Get(edg, r_1);   /* NOTE: this is always the same as cur_row, q_0 */
         lb_1  = bnd.lb;                     /* can't overflow the left edge */
         rb_1  = bnd.rb;                     /* can't overflow the right edge */

         for (t_0 = lb_1; t_0 < rb_1; t_0++) 
         {
            MMX3(qx1, t_0) = MY_Zero();
            IMX3(qx1, t_0) = MY_Zero();
            DMX3(qx1, t_0) = MY_Zero();
         }
      }

      /* SET CURRENT ROW TO PREVIOUS ROW */
      r_1b = r_0b;
      r_1e = r_0e;
   }

   /* UNROLLED FINAL ROW */
   /* At q_0 = 0, only N,B states are reachable. */
   q_0 = Q_range.beg;
   {
      q_1 = q_0 + 1;
      qx0 = q_0 % 2;
      qx1 = qx0 + 1;

      a = seq[q_0];
      A = AA_REV[a];

      /* get edgebound range */
      r_0b  = EDGEBOUNDS_GetIndex_byRow_Bck( edg, q_0 + 1 );
      r_0e  = EDGEBOUNDS_GetIndex_byRow_Bck( edg, q_0 );

      /* FINAL i = 0 row */
      a = seq[q_0];
      A = AA_REV[a];

      /* B STATE (SPARSE) */
      XMX(SP_B, q_0) = MY_Zero();
      for (r_1 = r_1b; r_1 > r_1e; r_1--) 
      {
         bnd   = EDGEBOUNDS_Get(edg, r_1);          /* bounds for current bound */
         lb_0  = MAX(bnd.lb, T_range.beg);          /* can't overflow the left edge */
         rb_0  = MIN(bnd.rb, T_range.end);          /* can't overflow the right edge */

         for (t_0 = rb_0-1; t_0 >= lb_0; t_0--)
         {
            t_1 = t_0 - 1;

            prv_sum  = XMX(SP_B, q_0);
            prv_M    = MY_Prod( MY_Prod(  MMX3(q_1, t_0), TSC(t_1, B2M) ), 
                                          MSC(t_0, A) );
            XMX(SP_B, q_0) = MY_Sum( prv_sum, prv_M );
         }
      }

      XMX(SP_J, q_0) = MY_Zero();
      XMX(SP_C, q_0) = MY_Zero();
      XMX(SP_E, q_0) = MY_Zero();

      prv_N = MY_Prod( XMX(SP_N, q_1), XSC(SP_N, SP_LOOP) );
      prv_B = MY_Prod( XMX(SP_B, q_0), XSC(SP_N, SP_MOVE) );
      XMX(SP_N, q_0) = MY_Sum( prv_N, prv_B );

      /* SCRUB FINAL ROW */
      for (r_1 = r_1b; r_1 > r_1e; r_1--) 
      {
         bnd   = EDGEBOUNDS_Get(edg, r_1);   /* NOTE: this is always the same as cur_row, q_0 */
         lb_1  = bnd.lb;                     /* can't overflow the left edge */
         rb_1  = bnd.rb;                     /* can't overflow the right edge */

         for (t_0 = lb_1; t_0 < rb_1; t_0++) 
         {
            MMX3(qx1, t_0) = MY_Zero();
            IMX3(qx1, t_0) = MY_Zero();
            DMX3(qx1, t_0) = MY_Zero();
         }
      }
   }
   
   #if DEBUG 
   {
      int dirty_vals = MATRIX_3D_Check_Clean( st_MX3 );
      printf("# POST-CLEAN CHECK -> BOUND BCK LIN:\t%d\n", dirty_vals);
      if (dirty_vals > 0) {
         MATRIX_3D_Dump( st_MX3, stdout );
      }
   }
   #endif

   sc_best     = XMX(SP_N, Q_range.beg);
   *sc_final   = sc_best;

   return STATUS_SUCCESS;
}

/* MATH RULES: These determine how probilities are summed and certain identities */

static 
inline 
float
MY_Sum(  const float    x,
         const float    y )
{
   return MATH_LogSum( x, y );
}

static 
inline
float
MY_Prod( const float    x,
         const float    y )
{
   return MATH_LogProd( x, y );
}

static 
inline
float
MY_Zero()
{
   return MATH_LogZero();
}

static 
inline
float 
MY_One()
{
   return MATH_LogOne();
}