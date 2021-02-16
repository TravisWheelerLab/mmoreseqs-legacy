/*******************************************************************************
 *     FILE:  bounded_fwdbck_linear.c
 *  PURPOSE:  Bounded Forward/Backward Algorithm 
 *            (Linear Space Alg)
 *
 *   AUTHOR:  Dave Rich
 *     BUGS:       
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
#include "../algs_linear/_algs_linear.h"
#include "../parsers/_parsers.h"

/* self header */
#include "_algs_sparse.h"
#include "bound_fwdbck_sparse.h"

/* NOTE: HOW TO CONVERT row-coords to diag-coords
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


/** FUNCTION:  run_Bound_Forward_Sparse()
 *  SYNOPSIS:  Perform Edge-Bounded Forward step of Cloud Search Algorithm.
 *             Runs traditional Forward-Backward Algorithm, but only performs
 *             computation on cells that fall within the bounds determined by
 *             the <edg> EDGEBOUNDS object, which stores a series of 
 *             (left-bound, right-bound) pairs sorted by row.  
 *             Normal state matrix is stored in linear space.
 *             <st_MX3> is size [3 * (Q + T + 1)]. Only requires size [2 * (T + 1)],
 *             but is reused from cloud_forward_().
 *             Final score produced by Forward is stored in <sc_final>.
 *
 *    RETURN:  Returns the final score of the Forward Algorithm.
 */
STATUS_FLAG 
run_Bound_Forward_Sparse_TEST(   const SEQUENCE*               query,         /* query sequence */
                                 const HMM_PROFILE*            target,        /* target HMM model */
                                 const int                     Q,             /* query length */
                                 const int                     T,             /* target length */
                                 MATRIX_3D_SPARSE* restrict    st_SMX_fwd,    /* normal state matrix */
                                 MATRIX_2D* restrict           sp_MX_fwd,     /* special state matrix */
                                 const EDGEBOUNDS*             edg,           /* edgebounds */
                                 const RANGE*                  dom_range,     /* (OPTIONAL) domain range for computing fwd/bck on specific domain. If NULL, computes complete fwd/bck. */
                                 float*                        sc_final )     /* (OUTPUT) final score */
{
   /* vars for matrix access for macros */
   MATRIX_3D_SPARSE*    st_SMX   = st_SMX_fwd;    /* normal state matrix */
   MATRIX_2D*           sp_MX    = sp_MX_fwd;     /* special state matrix */

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
   int      id_0;                            /* id in edgebound list (row/diag) */
   int      r_0;                             /* current index for current row */
   int      r_0b, r_0e;                      /* begin and end indices for current row in edgebound list */
   int      r_1;                             /* current index for previous row */
   int      r_1b, r_1e;                      /* begin and end indices for current row in edgebound list */
   int      le_0, re_0;                      /* right/left matrix bounds of current diag */
   int      lb_0, rb_0;                      /* bounds of current search space on current diag */
   int      lb_1, rb_1;                      /* bounds of current search space on previous diag */
   int      lb_2, rb_2;                      /* bounds of current search space on 2-back diag */
   int      lb_T, rb_T;                      /* truncated bounds within domain */

   /* vars for recurrance scores */
   float    prv_M, prv_I, prv_D;             /* previous (M) match, (I) insert, (D) delete states */
   float    prv_B, prv_E;                    /* previous (B) begin and (E) end states */
   float    prv_J, prv_N, prv_C;             /* previous (J) jump, (N) initial, and (C) terminal states */
   float    prv_loop, prv_move;              /* previous loop and move for special states */
   float    prv_sum, prv_best;               /* temp subtotaling vars */
   float    sc_best;                         /* final best scores */
   float    sc_M, sc_I, sc_D, sc_E;          /* match, insert, delete, end scores */

   /* vars for sparse matrix */
   MATRIX_3D_SPARSE*    mx;                     
   EDGEBOUNDS*          edg_inner;           /* edgebounds for search space of backward/forward */
   EDGEBOUNDS*          edg_outer;           /* edgebounds for sparse matrix shape */
   RANGE                T_range;             /* target range */
   RANGE                Q_range;             /* query range */
   bool                 is_q_0_in_dom_range; /* checks if current query position is inside the domain range */
   bool                 is_q_1_in_dom_range; /* checks if previous query position is inside the domain range */

   /* debugger tools */
   FILE*          dbfp;
   MATRIX_2D*     cloud_MX;
   MATRIX_2D*     cloud_MX3;
   MATRIX_3D*     test_MX;
   MATRIX_3D*     test_MX3;
   EDGEBOUNDS*    test_edg;
   int            num_writes;
   int            num_clears;

   /* initialize debugging matrix */
   #if DEBUG
   {
      cloud_MX    = debugger->cloud_MX;
      cloud_MX3   = debugger->cloud_MX3;
      test_MX     = debugger->test_MX;
      test_MX3    = debugger->test_MX3;
      test_edg    = debugger->test_edg;

      MATRIX_2D_Reuse( cloud_MX, Q+1, T+1 );
      MATRIX_2D_Fill( cloud_MX, 0.0 );
      int   num_cells   = MATRIX_2D_Cloud_Fill( cloud_MX, edg, 1.0 );
      int   num_cells2  = MATRIX_2D_Cloud_Count( cloud_MX );
      float pretotal    = MATRIX_2D_Total( cloud_MX );
      printf("Forward Pretotal = %d %d %f\n", num_cells, num_cells2, pretotal);

      // MATRIX_2D_Reuse( cloud_MX3, 3, (Q+1)+(T+1) );
      // MATRIX_2D_Fill( cloud_MX3, 0 );
      // MATRIX_3D_Reuse( test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
      // MATRIX_3D_Fill( test_MX, -INF );
      // MATRIX_3D_Reuse( test_MX3, NUM_NORMAL_STATES, 3, (Q+1)+(T+1) );
      // MATRIX_3D_Fill( test_MX3, -INF );
      // EDGEBOUNDS_Reuse( test_edg, Q, T );
      // MATRIX_2D_Fill( sp_MX, -INF );

      num_writes = 0;
      num_clears = 0;
   }
   #endif

   /* --------------------------------------------------------------------------------- */

   /* initialize logsum lookup table if it has not already been */
   MATH_Logsum_Init();
   // SEQUENCE_Digitize( query );

   /* query sequence */
   mx          = st_SMX;
   seq         = query->seq;
   N           = EDGEBOUNDS_GetSize( edg );
   /* local or global alignments? */
   is_local    = target->isLocal;
   sc_E        = (is_local) ? MY_One() : MY_Zero();
   if ( is_local == true ) {
      MY_One();
   }
   else {
      MY_Zero();
   }

   /* domain range (query sequence) */
   if (dom_range == NULL) {
      Q_range.beg = 0;
      Q_range.end = Q + 1;
   } else {
      Q_range = *dom_range;
   }
   /* target range */
   T_range.beg = 0;
   T_range.end = T + 1;

   /* UNROLLED INITIAL QUERY LOOP */
   q_0 = Q_range.beg;
   {   
      /* check if query position is in domain */
      is_q_0_in_dom_range = IS_IN_RANGE( Q_range.beg, Q_range.end, q_0 );
      /* get edgebound range */
      r_0b  = EDGEBOUNDS_GetIndex_byRow_Fwd( edg, q_0 );
      r_0e  = EDGEBOUNDS_GetIndex_byRow_Fwd( edg, q_0 + 1 );

      /* initialize special states */
      XMX(SP_E, q_0) = MY_Zero(); 
      XMX(SP_J, q_0) = MY_Zero(); 
      XMX(SP_C, q_0) = MY_Zero(); 
      /* S->N, p=1 */
      XMX(SP_N, q_0) = MY_One(); 
      /* S->N->B, no N-tail */
      XMX(SP_B, q_0) = XSC(SP_N,SP_MOVE); 

      /* only compute if in domain range */
      if ( is_q_0_in_dom_range == true )
      {
         /* FOR every BOUND in zero row */
         for (r_0 = r_0b; r_0 < r_0e; r_0++)
         {
            /* get bound data */
            bnd   = MATRIX_3D_SPARSE_GetBound_byIndex( st_SMX, r_0 );
            id_0  = bnd.id;
            lb_0  = MAX(bnd.lb - 1, T_range.beg);    /* can't overflow the left edge */
            rb_0  = MIN(bnd.rb, T_range.end);    /* can't overflow the right edge */

            /* fetch data mapping bound start location to data block in sparse matrix */
            qx0   = MATRIX_3D_SPARSE_GetOffset_ByIndex_Cur( st_SMX, r_0 );
            /* initial location for square matrix and mapping to sparse matrix */
            t_0   = lb_0;
            tx0   = t_0 - bnd.lb;    /* total_offset = offset_location - starting_location */

            float* restrict MTX_[NUM_NORMAL_STATES][1][1];
            MTX_[MAT_ST][0][0] = MATRIX_3D_SPARSE_GetX_byOffset( mx, qx0, tx0, MAT_ST );
            MTX_[INS_ST][0][0] = MATRIX_3D_SPARSE_GetX_byOffset( mx, qx0, tx0, INS_ST );
            MTX_[DEL_ST][0][0] = MATRIX_3D_SPARSE_GetX_byOffset( mx, qx0, tx0, DEL_ST );

            /* FOR every position in TARGET profile */
            for (t_0 = lb_0; t_0 < rb_0; t_0++, tx0++)
            {
               tx0   = t_0 - bnd.lb;

               *MTX_[M_ST][0][0]  = MY_Zero();
               *MTX_[I_ST][0][0]  = MY_Zero();
               *MTX_[D_ST][0][0]  = MY_Zero();

               /* zero column is -inf in logspace.  We can skip this step and convert to normal space now. */
               // MSMX(qx0, tx0) = MY_Zero();
               // ISMX(qx0, tx0) = MY_Zero();
               // DSMX(qx0, tx0) = MY_Zero(); 

               /* embed linear row into quadratic test matrix */
               #if DEBUG
               {
                  MX_2D(cloud_MX, q_0, t_0) += 2.0;
                  MX_3D(test_MX, MAT_ST, q_0, t_0) = MSMX(qx0, tx0);
                  MX_3D(test_MX, INS_ST, q_0, t_0) = ISMX(qx0, tx0);
                  MX_3D(test_MX, DEL_ST, q_0, t_0) = DSMX(qx0, tx0);
               }
               #endif

               MTX_[M_ST][0][0]   += NUM_NORMAL_STATES;
               MTX_[I_ST][0][0]   += NUM_NORMAL_STATES;
               MTX_[D_ST][0][0]   += NUM_NORMAL_STATES;
            }
         }
      }

      /* lookback one row */
      r_1b = r_0b;
      r_1e = r_0e;
   }
   
   /* MAIN QUERY LOOP */
   /* FOR every position in QUERY sequence (row in matrix) */
   for ( q_0 = Q_range.beg + 1; q_0 < Q_range.end; q_0++ )
   {
      q_1 = q_0 - 1;
      t_0 = 0;

      /* check if query position is in domain */
      is_q_0_in_dom_range = IS_IN_RANGE( Q_range.beg + 1, Q_range.end, q_0 );
      is_q_0_in_dom_range = (q_0 >= Q_range.beg && q_0 < Q_range.end);
      /* get edgebound range */
      r_0b  = EDGEBOUNDS_GetIndex_byRow_Fwd( edg, q_0 );
      r_0e  = EDGEBOUNDS_GetIndex_byRow_Fwd( edg, q_0 + 1 );

      /* Get next sequence character */
      a = seq[q_1];  /* off-by-one */
      A = AA_REV[a];
      // SEQUENCE_GetDigitAt( query, q_1 );

      /* Init E state for current row */
      XMX(SP_E, q_0) = MY_Zero();

      /* only compute if in domain range */
      if ( is_q_0_in_dom_range == true )
      {
         /* FOR every BOUND in current ROW */
         for (r_0 = r_0b; r_0 < r_0e; r_0++)
         {
            /* get bound data */
            bnd   = MATRIX_3D_SPARSE_GetBound_byIndex( st_SMX, r_0 );
            id_0  = bnd.id;
            lb_0  = MAX(bnd.lb - 1, T_range.beg);    /* can't overflow left edge. the leftmost cell will be set to -INF, so (-1) adds left padding cell.  */
            rb_0  = MIN(bnd.rb, T_range.end);        /* can't overflow right edge */

            /* fetch data mapping bound start location to data block in sparse matrix */
            qx0 = MATRIX_3D_SPARSE_GetOffset_ByIndex_Cur( st_SMX, r_0 );
            qx1 = MATRIX_3D_SPARSE_GetOffset_ByIndex_Prv( st_SMX, r_0 );
            /* initial location for square matrix and mapping to sparse matrix */
            t_0 = lb_0;
            tx0 = t_0 - bnd.lb;    /* total_offset = offset_location - starting_location */
            tx1 = tx0 - 1;

            float* restrict MTX_[NUM_NORMAL_STATES][2][2];
            int st_[] = { MAT_ST, INS_ST, DEL_ST };
            int qx_[] = { qx0, qx1 };
            int tx_[] = { tx0, tx1 };
            for (int i = 0; i < 2; i++) {
               for (int j = 0; j < 2; j++) {
                  for (int k = 0; k < NUM_NORMAL_STATES; k++) {
                     MTX_[k][i][j] = MATRIX_3D_SPARSE_GetX_byOffset( mx, qx_[i], tx_[j], st_[k] );
                  }
               }
            }

            /* UNROLLED INITIAL TARGET LOOP: special case for left edge of range */
            t_0 = lb_0;
            {
               tx0 = t_0 - bnd.lb;

               *MTX_[M_ST][0][0]    += MY_Zero();
               *MTX_[I_ST][0][0]    += MY_Zero();
               *MTX_[D_ST][0][0]    += MY_Zero();

               #if FALSE
               {
                  MSMX(qx0, tx0) = MY_Zero();
                  ISMX(qx0, tx0) = MY_Zero();
                  DSMX(qx0, tx0) = MY_Zero();
               }
               #endif

               /* embed linear row into quadratic test matrix */
               #if DEBUG
               {
                  MX_2D(cloud_MX, q_0, t_0) += 2.0;
                  MX_3D(test_MX, MAT_ST, q_0, t_0) = MSMX(qx0, tx0);
                  MX_3D(test_MX, INS_ST, q_0, t_0) = ISMX(qx0, tx0);
                  MX_3D(test_MX, DEL_ST, q_0, t_0) = DSMX(qx0, tx0);
               }
               #endif

               for (int i = 0; i < 2; i++) {
                  for (int j = 0; j < 2; j++) {
                     for (int k = 0; k < NUM_NORMAL_STATES; k++) {
                        MTX_[k][i][j] += NUM_NORMAL_STATES;
                     }
                  }
               }
            }

            

            /* MAIN TARGET LOOP */
            /* FOR every position in TARGET profile */
            for ( t_0 = lb_0 + 1; t_0 < rb_0 - 1; t_0++ )
            {
               t_1 = t_0 - 1; 
               tx0 = t_0 - bnd.lb;
               tx1 = tx0 - 1;

               /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
               /* best previous state transition (match takes the diag element of each prev state) */
               prv_M    = MY_Prod( *MTX_[M_ST][1][1], TSC(t_1, M2M) );
               prv_I    = MY_Prod( *MTX_[I_ST][1][1], TSC(t_1, I2M) );
               prv_D    = MY_Prod( *MTX_[D_ST][1][1], TSC(t_1, TM) );
               prv_B    = MY_Prod( XMX(SP_B, q_1), TSC(t_1, B2M) ); /* from begin match state (new alignment) */
               /* best-to-match */
               prv_sum  = MY_Sum( MY_Sum( prv_M, prv_I ),
                                  MY_Sum( prv_B, prv_D ) );
               *MTX_[M_ST][0][0] = MY_Prod( prv_sum, MSC(t_0, A) );

               /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
               /* previous states (match takes the previous row (upper) of each state) */
               prv_M    = MY_Prod( *MTX_[M_ST][1][0], TSC(t_0, M2I) );
               prv_I    = MY_Prod( *MTX_[I_ST][1][0], TSC(t_0, I2I) );
               /* best-to-insert */
               prv_sum  = MY_Sum( prv_M, prv_I );
               *MTX_[I_ST][0][0] = MY_Prod( prv_sum, ISC(t_0, A) );

               /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) */
               /* previous states (match takes the previous column (left) of each state) */
               prv_M    = MY_Prod( *MTX_[M_ST][0][1], TSC(t_1, M2D) );
               prv_D    = MY_Prod( *MTX_[D_ST][0][1], TSC(t_1, TD) );
               /* best-to-delete */
               prv_sum  = MY_Sum( prv_M, prv_D );
               *MTX_[D_ST][0][0] = prv_sum;

               /* UPDATE E STATE */
               prv_M    = MY_Prod( *MTX_[M_ST][0][0], sc_E );
               prv_D    = MY_Prod( *MTX_[D_ST][0][0], sc_E );
               prv_E    = XMX(SP_E, q_0);
               /* best-to-e-state */
               prv_sum  = MY_Sum( MY_Sum( prv_M, prv_D ),
                                                 prv_E );
               XMX(SP_E, q_0) = prv_sum;

               #if FALSE
               {
                  /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
                  /* best previous state transition (match takes the diag element of each prev state) */
                  prv_M    = MY_Prod( MSMX(qx1, tx1), TSC(t_1, M2M) );
                  prv_I    = MY_Prod( ISMX(qx1, tx1), TSC(t_1, I2M) );
                  prv_D    = MY_Prod( DSMX(qx1, tx1), TSC(t_1, TM) );
                  prv_B    = MY_Prod( XMX(SP_B, q_1), TSC(t_1, B2M) ); /* from begin match state (new alignment) */
                  /* best-to-match */
                  prv_sum  = MY_Sum( MY_Sum( prv_M, prv_I ),
                                    MY_Sum( prv_B, prv_D ) );
                  MSMX(qx0, tx0) = MY_Prod( prv_sum, MSC(t_0, A) );

                  /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
                  /* previous states (match takes the previous row (upper) of each state) */
                  prv_M    = MY_Prod( MSMX(qx1, tx0), TSC(t_0, M2I) );
                  prv_I    = MY_Prod( ISMX(qx1, tx0), TSC(t_0, I2I) );
                  /* best-to-insert */
                  prv_sum  = MY_Sum( prv_M, prv_I );
                  ISMX(qx0, tx0) = MY_Prod( prv_sum, ISC(t_0, A) );

                  /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) */
                  /* previous states (match takes the previous column (left) of each state) */
                  prv_M    = MY_Prod( MSMX(qx0, tx1), TSC(t_1, M2D) );
                  prv_D    = MY_Prod( DSMX(qx0, tx1), TSC(t_1, TD) );
                  /* best-to-delete */
                  prv_sum  = MY_Sum( prv_M, prv_D );
                  DSMX(qx0, tx0) = prv_sum;

                  /* UPDATE E STATE */
                  prv_M    = MY_Prod( MSMX(qx0, tx0), sc_E );
                  prv_D    = MY_Prod( DSMX(qx0, tx0), sc_E );
                  prv_E    = XMX(SP_E, q_0);
                  /* best-to-e-state */
                  prv_sum  = MY_Sum( MY_Sum( prv_M, prv_D ),
                                                   prv_E );
                  XMX(SP_E, q_0) = prv_sum;
               }
               #endif

               /* embed linear row into quadratic test matrix */
               #if DEBUG
               {
                  MX_2D(cloud_MX, q_0, t_0) += 2.0;
                  MX_3D(test_MX, MAT_ST, q_0, t_0) = MSMX(qx0, tx0);
                  MX_3D(test_MX, INS_ST, q_0, t_0) = ISMX(qx0, tx0);
                  MX_3D(test_MX, DEL_ST, q_0, t_0) = DSMX(qx0, tx0);
               }
               #endif

               for (int i = 0; i < 2; i++) {
                  for (int j = 0; j < 2; j++) {
                     for (int k = 0; k < NUM_NORMAL_STATES; k++) {
                        MTX_[k][i][j] += NUM_NORMAL_STATES;
                     }
                  }
               }
            }

            /* UNROLLED FINAL TARGET LOOP: special case for right edge of range (only when range is greater than one cell) */
            if ( t_0 = rb_0 - 1, (rb_0 - lb_0) > 1 )  
            {
               t_1 = t_0 - 1; 
               tx0 = t_0 - bnd.lb;
               tx1 = tx0 - 1;

               /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
               /* best previous state transition (match takes the diag element of each prev state) */
               prv_M    = MY_Prod( MSMX(qx1, tx1), TSC(t_1, M2M) );
               prv_I    = MY_Prod( ISMX(qx1, tx1), TSC(t_1, I2M) );
               prv_D    = MY_Prod( DSMX(qx1, tx1), TSC(t_1, TM) );
               prv_B    = MY_Prod( XMX(SP_B, q_1), TSC(t_1, B2M) );    /* from begin match state (new alignment) */
               /* sum-to-match */
               prv_sum  = MY_Sum( MY_Sum( prv_M, prv_I ),
                                  MY_Sum( prv_D, prv_B ) );
               MSMX(qx0, tx0) = MY_Prod( prv_sum, MSC(t_0, A) );

               /* FIND SUM OF PATHS TO INSERT STATE (unrolled) */
               ISMX(qx0, tx0) = MY_Zero();

               /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) (unrolled) */
               /* previous states (match takes the left element of each state) */
               prv_M    = MY_Prod( MSMX(qx0, tx1), TSC(t_1, M2D) );
               prv_D    = MY_Prod( DSMX(qx0, tx1), TSC(t_1, TD) );
               /* sum-to-delete */
               prv_sum  = MY_Sum( prv_M, prv_D );
               DSMX(qx0, tx0) = prv_sum;

               /* UPDATE E STATE (unrolled) */
               prv_E    = XMX(SP_E, q_0);
               prv_M    = MSMX(qx0, tx0);
               prv_D    = DSMX(qx0, tx0);
               /* best-to-begin */
               prv_sum  = MY_Sum( MY_Sum( prv_D, prv_M ),
                                          prv_E ); 
               XMX(SP_E, q_0) = prv_sum;

               /* embed linear row into quadratic test matrix */
               #if DEBUG
               {
                  MX_2D(cloud_MX, q_0, t_0) += 2.0;
                  MX_3D(test_MX, MAT_ST, q_0, t_0) = MSMX(qx0, tx0);
                  MX_3D(test_MX, INS_ST, q_0, t_0) = ISMX(qx0, tx0);
                  MX_3D(test_MX, DEL_ST, q_0, t_0) = DSMX(qx0, tx0);
               }
               #endif
            }
         }
      }
      
      /* SPECIAL STATES */
      /* J state */
      prv_J    = MY_Prod( XMX(SP_J, q_1), XSC(SP_J, SP_LOOP) );       /* J->J */
      prv_E    = MY_Prod( XMX(SP_E, q_0), XSC(SP_E, SP_LOOP) );       /* E->J is E's "loop" */
      prv_sum  = MY_Sum( prv_J, prv_E );
      XMX(SP_J, q_0) = prv_sum;

      /* C state */
      prv_C    = MY_Prod( XMX(SP_C, q_1), XSC(SP_C, SP_LOOP) );
      prv_E    = MY_Prod( XMX(SP_E, q_0), XSC(SP_E, SP_MOVE) );
      prv_sum  = MY_Sum( prv_C, prv_E );
      XMX(SP_C, q_0) = prv_sum;

      /* N state */
      prv_N    = MY_Prod( XMX(SP_N, q_1), XSC(SP_N, SP_LOOP) );
      XMX(SP_N, q_0) = prv_N;

      /* B state */
      prv_N    = MY_Prod( XMX(SP_N, q_0), XSC(SP_N, SP_MOVE) );         /* N->B is N's move */
      prv_J    = MY_Prod( XMX(SP_J, q_0), XSC(SP_J, SP_MOVE) );         /* J->B is J's move */
      prv_sum  = MY_Sum( prv_N, prv_J ); 
      XMX(SP_B, q_0) = prv_sum;

      /* SET CURRENT ROW TO PREVIOUS ROW */
      r_1b = r_0b;
      r_1e = r_0e;
   }

   /* T state */
   sc_best     = MY_Prod( XMX(SP_C, Q_range.end - 1), XSC(SP_C, SP_MOVE) );
   *sc_final   = sc_best;

   /* output test matrix */
   #if DEBUG
   {
      // int   num_cells2  = MATRIX_2D_Cloud_Count( cloud_MX );
      // float pretotal    = MATRIX_2D_Total( cloud_MX );
      // printf("Forward Pretotal = %d %f\n", num_cells2, pretotal);

      // printf("POST => travis testing...\n");
      // int truetot = (Q+1) * (T+1);
      // int subtot = 0;
      // float numtotal = 0.0f;
      // int hit, non, missed, doublehit, badhit;
      // for (int i = 0; i < cloud_MX->R; i++) {
      //    for (int j = 0; j < cloud_MX->C; j++) {
      //       float val = MX_2D( cloud_MX, i, j );
      //       numtotal += val;
      //       if ( val == 0.0 ) non += 1;
      //       elif ( val == 1.0 ) missed += 1;
      //       elif ( val == 3.0 ) hit += 1;
      //       elif ( (int)val % 2 == 0 ) badhit += 1;
      //       else {
      //          doublehit += 1;
      //          printf("DOUBLE HIT AT: (%d,%d)\n", i, j);
      //       } 
      //    }
      // }
      // printf("NUMTOTAL: %f\n", numtotal);
      // subtot = non + missed + hit + doublehit;
      // printf("TRAVIS_COUNTS: hit= %d, non= %d, missed= %d, double= %d, bad= %d, totals= %d ? %d\n", 
      //    hit, non, missed, doublehit, badhit, subtot, truetot );
      // MATRIX_2D_Save( cloud_MX, "test_output/traviscount_fwd.mx");
   }
   #endif

   return STATUS_SUCCESS;
}


/** FUNCTION:  run_Bound_Forward_Sparse()
 *  SYNOPSIS:  Perform Edge-Bounded Forward step of Cloud Search Algorithm.
 *             Runs traditional Forward-Backward Algorithm, but only performs
 *             computation on cells that fall within the bounds determined by
 *             the <edg> EDGEBOUNDS object, which stores a series of 
 *             (left-bound, right-bound) pairs sorted by row.  
 *             Normal state matrix is stored in linear space.
 *             <st_MX3> is size [3 * (Q + T + 1)]. Only requires size [2 * (T + 1)],
 *             but is reused from cloud_forward_().
 *             Final score produced by Forward is stored in <sc_final>.
 *
 *    RETURN:  Returns the final score of the Forward Algorithm.
 */
STATUS_FLAG 
run_Bound_Backward_Sparse_TEST(  const SEQUENCE*               query,         /* query sequence */
                                 const HMM_PROFILE*            target,        /* target HMM model */
                                 const int                     Q,             /* query length */
                                 const int                     T,             /* target length */
                                 MATRIX_3D_SPARSE* restrict    st_SMX_bck,    /* normal state matrix */
                                 MATRIX_2D* restrict           sp_MX_bck,     /* special state matrix */
                                 const EDGEBOUNDS*             edg,           /* edgebounds */
                                 const RANGE*                  dom_range,     /* (OPTIONAL) domain range for computing fwd/bck on specific domain. If NULL, computes complete fwd/bck. */
                                 float*                        sc_final )     /* (OUTPUT) final score */
{
   /* vars for matrix access for macros */
   MATRIX_3D_SPARSE*    st_SMX   = st_SMX_bck;        /* normal state matrix */
   MATRIX_2D*           sp_MX    = sp_MX_bck;         /* special state matrix */

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
   int      lb_T, rb_T;                      /* truncated bounds within domain */

   /* vars for recurrance scores */
   float    prv_M, prv_I, prv_D;             /* previous (M) match, (I) insert, (D) delete states */
   float    prv_B, prv_E;                    /* previous (B) begin and (E) end states */
   float    prv_J, prv_N, prv_C;             /* previous (J) jump, (N) initial, and (C) terminal states */
   float    prv_loop, prv_move;              /* previous loop and move for special states */
   float    prv_sum, prv_best;               /* temp subtotaling vars */
   float    sc_best;                         /* final best scores */
   float    sc_M, sc_I, sc_D, sc_E;          /* match, insert, delete, end scores */

   /* vars for sparse matrix */
   MATRIX_3D_SPARSE*   mx;
   EDGEBOUNDS*          edg_inner;           /* edgebounds for search space of backward/forward */
   EDGEBOUNDS*          edg_outer;           /* edgebounds for sparse matrix shape */
   RANGE                T_range;             /* target range */
   RANGE                Q_range;             /* query range */
   bool                 is_q_0_in_dom_range; /* checks if current query position is inside the domain range */
   bool                 is_q_1_in_dom_range; /* checks if previous query position is inside the domain range */

   /* debugger tools */
   FILE*                dbfp;
   MATRIX_2D*           cloud_MX;
   MATRIX_2D*           cloud_MX3;
   MATRIX_3D*           test_MX;
   MATRIX_3D*           test_MX3;
   EDGEBOUNDS*          test_edg;
   int                  num_writes;
   int                  num_clears;

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

   /* --------------------------------------------------------------------------------- */

   /* initialize logsum lookup table if it has not already been */
   MATH_Logsum_Init();
   // SEQUENCE_Digitize( query );

   /* query sequence */
   mx          = st_SMX;
   seq         = query->seq;
   N           = EDGEBOUNDS_GetSize( edg );
   /* local or global alignments? */
   is_local    = target->isLocal;
   sc_E        = (is_local) ? MY_One() : MY_Zero();
   if ( is_local == true ) {
      MY_One();
   }
   else {
      MY_Zero();
   }

   /* domain range (query sequence) */
   if (dom_range == NULL) {
      Q_range.beg = 0;
      Q_range.end = Q;
   } else {
      Q_range = *dom_range;
   }
   /* valid target range */
   T_range.beg = 1;
   T_range.end = T;

   /* UNROLLED INITIAL QUERY LOOP */
   q_0 = Q_range.end;
   {
      /* if inside domain */
      is_q_0_in_dom_range = IS_IN_RANGE( Q_range.beg, Q_range.end, q_0 );
      /* get edgebound range */
      r_0b  = EDGEBOUNDS_GetIndex_byRow_Bck( edg, q_0 + 1 );
      r_0e  = EDGEBOUNDS_GetIndex_byRow_Bck( edg, q_0 );

      /* INIT SPECIAL STATES */
      XMX(SP_J, q_0) = MY_Zero();
      XMX(SP_B, q_0) = MY_Zero();
      XMX(SP_N, q_0) = MY_Zero();
      XMX(SP_C, q_0) = XSC(SP_C, SP_MOVE);
      XMX(SP_E, q_0) = MY_Prod( XMX(SP_C, q_0), XSC(SP_E, SP_MOVE) );

      /* if sequence position is in domain range */ 
      // if ( is_q_0_in_dom_range == true )
      {
         bnd   = MATRIX_3D_SPARSE_GetBound_byIndex( st_SMX, r_0b );

         /* FOR every SPAN in current ROW */
         for (r_0 = r_0b; r_0 > r_0e; r_0--) 
         {
            /* get bound data */
            // bnd   = EDG_X(edg, r_0);            /* bounds for current bound */
            bnd   = MATRIX_3D_SPARSE_GetBound_byIndex( st_SMX, r_0 );
            lb_0  = MAX(bnd.lb, T_range.beg);      /* can't overflow left edge */
            rb_0  = MIN(bnd.rb, T_range.end);  /* can't overflow right edge */

            /* fetch data mapping bound start location to data block in sparse matrix */
            qx0   = MATRIX_3D_SPARSE_GetOffset_ByIndex_Cur( st_SMX, r_0 );
            qx1   = MATRIX_3D_SPARSE_GetOffset_ByIndex_Nxt( st_SMX, r_0 );
            /* location for square matrix and mapping to sparse matrix */
            t_0 = lb_0;
            tx0 = t_0 - bnd.lb;

            /* UNROLLED INITIAL TARGET LOOP */
            t_0 = rb_0;
            {
               t_1 = t_0 + 1;
               tx0 = t_0 - bnd.lb;
               tx1 = tx0 + 1;

               MSMX(qx0, tx0) = XMX(SP_E, q_0);
               ISMX(qx0, tx0) = MY_Zero();
               DSMX(qx0, tx0) = XMX(SP_E, q_0);

               #if DEBUG 
               {
                  MX_2D(cloud_MX, q_0, t_0) += 2.0;
                  MX_3D(test_MX, MAT_ST, q_0, t_0) = MSMX(qx0, tx0);
                  MX_3D(test_MX, INS_ST, q_0, t_0) = ISMX(qx0, tx0);
                  MX_3D(test_MX, DEL_ST, q_0, t_0) = DSMX(qx0, tx0);
               }
               #endif
            }

            /* MAIN TARGET LOOP */
            /* FOR every position of TARGET in SPAN */
            for (t_0 = rb_0 - 1; t_0 >= lb_0; t_0--)
            {
               t_1 = t_0 + 1;
               tx0 = t_0 - bnd.lb;
               tx1 = tx0 + 1;

               prv_E = MY_Prod( XMX(SP_E, Q), sc_E );
               prv_D = MY_Prod( DSMX(qx0, tx1), TSC(t_0, M2D) );
               MSMX(qx0, tx0) = MY_Sum( prv_E, prv_D );

               ISMX(qx0, tx0) = MY_Zero();

               prv_E = MY_Prod( XMX(SP_E, Q), sc_E );
               prv_D = MY_Prod( DSMX(qx0, tx1), TSC(t_0, TD) );
               DSMX(qx0, tx0) = MY_Sum( prv_E, prv_D );

               #if DEBUG 
               {
                  MX_2D(cloud_MX, q_0, t_0) += 2.0;
                  MX_3D(test_MX, MAT_ST, q_0, t_0) = MSMX(qx0, tx0);
                  MX_3D(test_MX, INS_ST, q_0, t_0) = ISMX(qx0, tx0);
                  MX_3D(test_MX, DEL_ST, q_0, t_0) = DSMX(qx0, tx0);
               }  
               #endif
            }
         }
      }

      /* init lookback 1 row */
      r_1b = r_0b;
      r_1e = r_0e;
   }
   
   /* MAIN QUERY LOOP */
   /* FOR every position in QUERY */
   for ( q_0 = Q_range.end - 1; q_0 > Q_range.beg; q_0-- )
   {
      q_1 = q_0 + 1;

      /* if inside domain */
      is_q_0_in_dom_range = IS_IN_RANGE( Q_range.beg, Q_range.end, q_0 );
      is_q_1_in_dom_range = IS_IN_RANGE( Q_range.beg, Q_range.end, q_1 );
      /* get edgebound range */
      r_0b  = EDGEBOUNDS_GetIndex_byRow_Bck( edg, q_0 + 1 );
      r_0e  = EDGEBOUNDS_GetIndex_byRow_Bck( edg, q_0 );

      /* Get next sequence character */
      a = seq[q_0];
      A = AA_REV[a];
      // SEQUENCE_GetDigitAt( query, q_0 );

      /* UPDATE B STATE */
      XMX(SP_B, q_0) = MY_Zero();
      /* if previous q is in domain range, update B state */
      // if ( is_q_1_in_dom_range == true )
      {
         for (r_1 = r_1b; r_1 > r_1e; r_1--) 
         {
            /* get bound data */
            // bnd   = EDG_X(edg, r_1);              /* bounds for current bound */
            bnd   = MATRIX_3D_SPARSE_GetBound_byIndex( st_SMX, r_1 );
            lb_0  = MAX(bnd.lb, T_range.beg);   /* can't overflow left edge */
            rb_0  = MIN(bnd.rb, T_range.end);   /* can't overflow right edge */

            /* fetch data location to bound start location (in offset) */
            // qx1 = VECTOR_INT_Get( st_SMX->imap_cur, r_1 );    /* (q_0, t_0) location offset */
            qx1   = MATRIX_3D_SPARSE_GetOffset_ByIndex_Cur( st_SMX, r_1 );
            /* location for square matrix and mapping to sparse matrix */
            t_0 = lb_0;
            tx0 = t_0 - bnd.lb;    /* total_offset = offset_location - starting_location */

            for (t_0 = rb_0 - 1; t_0 >= lb_0; t_0--)
            {
               t_1 = t_0 - 1;
               tx0 = t_0 - bnd.lb;
               tx1 = tx0 - 1;

               prv_sum  = XMX(SP_B, q_0);
               prv_M    = MY_Prod( MY_Prod( MSMX(qx1, tx0), TSC(t_1, B2M) ), 
                                            MSC(t_0, A) );
               XMX(SP_B, q_0) = MY_Sum( prv_sum, prv_M);
            }
         }
      }

      /* UPDATE SPECIAL STATES */
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

      // if ( is_q_0_in_dom_range == true )
      {
         /* FOR every SPAN in current ROW */
         for (r_0 = r_0b; r_0 > r_0e; r_0--)
         {
            /* get bound data */
            // bnd   = EDG_X(edg, r_0);       /* bounds for current bound */
            bnd   = MATRIX_3D_SPARSE_GetBound_byIndex( st_SMX, r_0 );
            lb_0  = MAX(bnd.lb, T_range.beg);         /* can't overflow left edge */
            rb_0  = MIN(bnd.rb, T_range.end);     /* can't overflow right edge */
            
            /* fetch data location to bound start location (in offset) */
            // qx0   = VECTOR_INT_Get( st_SMX->imap_cur, r_0 );    /* (q_0, t_0) location offset */
            // qx1   = VECTOR_INT_Get( st_SMX->imap_nxt, r_0 );    /* (q_0, t_0) location offset */
            qx0   = MATRIX_3D_SPARSE_GetOffset_ByIndex_Cur( st_SMX, r_0 );
            qx1   = MATRIX_3D_SPARSE_GetOffset_ByIndex_Nxt( st_SMX, r_0 );
            /* location for square matrix and mapping to sparse matrix */
            t_0 = T;
            tx0 = t_0 - bnd.lb;    /* total_offset = offset_location - starting_location */

            /* UNROLLED INITIAL TARGET LOOP */
            t_0 = rb_0;
            {
               t_1 = t_0 + 1;
               tx0 = t_0 - bnd.lb;
               tx1 = tx0 + 1;

               MSMX(qx0, tx0) = XMX(SP_E, q_0);
               ISMX(qx0, tx0) = MY_Zero();
               DSMX(qx0, tx0) = XMX(SP_E, q_0);

               #if DEBUG 
               {
                  MX_2D(cloud_MX, q_0, t_0) += 2.0;
                  MX_3D(test_MX, MAT_ST, q_0, t_0) = MSMX(qx0, tx0);
                  MX_3D(test_MX, INS_ST, q_0, t_0) = ISMX(qx0, tx0);
                  MX_3D(test_MX, DEL_ST, q_0, t_0) = DSMX(qx0, tx0);
               }
               #endif
            }

            /* MAIN TARGET LOOP */
            /* FOR every position of TARGET in SPAN */
            for (t_0 = rb_0 - 1; t_0 >= lb_0; t_0--)
            {
               t_1 = t_0 + 1;
               tx0 = t_0 - bnd.lb;
               tx1 = tx0 + 1;

               /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
               prv_M    = MY_Prod( MY_Prod(  MSMX(qx1, tx1), TSC(t_0, M2M) ), 
                                             MSC(t_1, A) );
               prv_I    = MY_Prod( MY_Prod(  ISMX(qx1, tx0), TSC(t_0, M2I) ), 
                                             ISC(t_1, A) );
               prv_D    = MY_Prod( DSMX(qx0, tx1), TSC(t_0, M2D) );
               prv_E    = MY_Prod( XMX(SP_E, q_0), sc_E );     /* from end match state (new alignment) */
               /* best-to-match */
               prv_sum  = MY_Sum( MY_Sum( prv_M, prv_I ),
                                  MY_Sum( prv_E, prv_D ) );
               MSMX(qx0, tx0) = prv_sum;

               /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
               prv_M    = MY_Prod( MY_Prod( MSMX(qx1, tx1), TSC(t_0, I2M) ),
                                            MSC(t_1, A) );
               prv_I    = MY_Prod( MY_Prod( ISMX(qx1, tx0), TSC(t_0, I2I) ), 
                                            ISC(t_0, A) );
               /* best-to-insert */
               prv_sum  = MY_Sum( prv_M, prv_I );
               ISMX(qx0, tx0) = prv_sum;

               /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
               prv_M    = MY_Prod( MY_Prod( MSMX(qx1, tx1), TSC(t_0, TM) ), 
                                            MSC(t_1, A) );
               prv_D    = MY_Prod( DSMX(qx0, tx1), TSC(t_0, TD) );
               prv_E    = MY_Prod( XMX(SP_E, q_0), sc_E );
               /* best-to-delete */
               prv_sum  = MY_Sum( prv_M, 
                          MY_Sum( prv_D, prv_E ) );
               DSMX(qx0, tx0) = prv_sum;

               #if DEBUG 
               {
                  MX_2D(cloud_MX, q_0, t_0) += 2.0;
                  MX_3D(test_MX, MAT_ST, q_0, t_0) = MSMX(qx0, tx0);
                  MX_3D(test_MX, INS_ST, q_0, t_0) = ISMX(qx0, tx0);
                  MX_3D(test_MX, DEL_ST, q_0, t_0) = DSMX(qx0, tx0);
               }
               #endif
            }
         }
      }
      
      /* SET CURRENT ROW TO PREVIOUS ROW */
      r_1b = r_0b;
      r_1e = r_0e;
   }

   /* UNROLLED FINAL ROW */
   q_0 = Q_range.beg;
   {
      /* At q_0 = 0, only N,B states are reachable. */
      q_1 = q_0 + 1;

      /* get edgebound range */
      r_0b  = EDGEBOUNDS_GetIndex_byRow_Bck( edg, q_0 + 1 );
      r_0e  = EDGEBOUNDS_GetIndex_byRow_Bck( edg, q_0 );

      /* FINAL i = 0 row */
      a = seq[q_0];
      A = AA_REV[a];
      // SEQUENCE_GetDigitAt( query, q_0 );

      /* UPDATE B STATE */
      XMX(SP_B, q_0) = MY_Zero();
      /* if previous q is in domain, update B state */
      // if ( is_q_1_in_dom_range == true )
      {
         for (r_1 = r_1b; r_1 > r_1e; r_1--) 
         {
            // bnd   = EDG_X(edg, r_1);       /* bounds for current bound */
            bnd   = MATRIX_3D_SPARSE_GetBound_byIndex( st_SMX, r_1 );
            lb_0  = MAX(bnd.lb, T_range.beg);   /* can't overflow left edge */
            rb_0  = MIN(bnd.rb, T_range.end);   /* can't overflow right edge */

            /* fetch data location to bound start location (in offset) */
            // qx1 = VECTOR_INT_Get( st_SMX->imap_cur, r_1 );    /* (q_0, t_0) location offset */
            qx1 = MATRIX_3D_SPARSE_GetOffset_ByIndex_Cur( st_SMX, r_1 );

            /* location for square matrix and mapping to sparse matrix */
            t_0 = lb_0;
            tx0 = t_0 - bnd.lb;    /* total_offset = offset_location - starting_location */

            for (t_0 = rb_0 - 1; t_0 >= lb_0; t_0--)
            {
               t_1 = t_0 - 1;
               /* calculate offset from beginning of sparse data block */ 
               tx0 = t_0 - bnd.lb;
               tx1 = tx0 - 1;

               prv_sum  = XMX(SP_B, q_0);
               prv_M    = MY_Prod( MY_Prod(  MSMX(qx1, tx0), TSC(t_1, B2M) ), 
                                             MSC(t_0, A) );
               XMX(SP_B, q_0) = MY_Sum( prv_sum, prv_M );
            }
         }
      }

      /* UPDATE SPECIAL STATES */
      XMX(SP_J, q_0) = MY_Zero();
      XMX(SP_C, q_0) = MY_Zero();
      XMX(SP_E, q_0) = MY_Zero();

      prv_N = MY_Prod( XMX(SP_N, q_1), XSC(SP_N, SP_LOOP) );
      prv_B = MY_Prod( XMX(SP_B, q_0), XSC(SP_N, SP_MOVE) );
      XMX(SP_N, q_0) = MY_Sum( prv_N, prv_B );
   }

   sc_best     = XMX(SP_N, Q_range.beg);
   *sc_final   = sc_best;

   #if DEBUG
   {

   }
   #endif

   return STATUS_SUCCESS;
}

/* MATH RULES: These determine how probilities are summed and certain identities */

static 
inline 
float
MY_Sum(  const float    x,
         const float    y )
{
   return MATH_NormalSum( x, y );
}

static 
inline
float
MY_Prod( const float    x,
         const float    y )
{
   return MATH_NormalProd( x, y );
}

static 
inline
float
MY_Zero()
{
   return MATH_NormalZero();
}

static 
inline
float 
MY_One()
{
   return MATH_NormalOne();
}
