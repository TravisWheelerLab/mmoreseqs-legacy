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
#include "../utilities/utilities.h"
#include "../objects/objects.h"
#include "../algs_linear/algs_linear.h"
#include "../parsers/parsers.h"

/* self header */
#include "_algs_sparse.h"
#include "bound_fwdbck_sparse.h"

/* NOTE: HOW TO CONVERT row-coords to diag-coords
 *       MMX3(i-1,j-1) => MMX3(, d_2)
 *       MMX3(i,  j-1) => MMX3(, d_1)
 *       MMX3(i,  j  ) => MMX3(, d_1)
 */

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
int run_Bound_Forward_Sparse(    const SEQUENCE*      query,         /* query sequence */
                                 const HMM_PROFILE*   target,        /* target HMM model */
                                 const int            Q,             /* query length */
                                 const int            T,             /* target length */
                                 MATRIX_3D_SPARSE*    st_SMX,         /* normal state matrix */
                                 MATRIX_2D*           sp_MX,         /* special state matrix */
                                 EDGEBOUNDS*          edg,           /* edgebounds */
                                 RANGE*               dom_range,     /* (OPTIONAL) domain range for computing fwd/bck on specific domain. If NULL, computes complete fwd/bck. */
                                 float*               sc_final )     /* (OUTPUT) final score */
{
   /* vars for matrix access */
   MATRIX_3D_SPARSE*    n_mx;                /* normal state matrix */
   MATRIX_2D*           s_mx;                /* special state matrix */

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
   BOUND*   bnd;                             /* current bound */
   int      id;                              /* id in edgebound list (row/diag) */
   int      r_0;                             /* current index for current row */
   int      r_0b, r_0e;                      /* begin and end indices for current row in edgebound list */
   int      r_1;                             /* current index for previous row */
   int      r_1b, r_1e;                      /* begin and end indices for current row in edgebound list */
   int      le_0, re_0;                      /* right/left matrix bounds of current diag */
   int      lb_0, rb_0;                      /* bounds of current search space on current diag */
   int      lb_1, rb_1;                      /* bounds of current search space on previous diag */
   int      lb_2, rb_2;                      /* bounds of current search space on 2-back diag */
   bool     lb_T, rb_T;                      /* checks if touches edge of matrix */

   /* vars for recurrance scores */
   float    prv_M, prv_I, prv_D;             /* previous (M) match, (I) insert, (D) delete states */
   float    prv_B, prv_E;                    /* previous (B) begin and (E) end states */
   float    prv_J, prv_N, prv_C;             /* previous (J) jump, (N) initial, and (C) terminal states */
   float    prv_loop, prv_move;            /* previous loop and move for special states */
   float    prv_sum, prv_best;             /* temp subtotaling vars */
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
      MATRIX_2D_Fill( cloud_MX, 0 );
      MATRIX_2D_Reuse( cloud_MX3, 3, (Q+1)+(T+1) );
      MATRIX_2D_Fill( cloud_MX3, 0 );
      MATRIX_3D_Reuse( test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
      MATRIX_3D_Fill( test_MX, -INF );
      MATRIX_3D_Reuse( test_MX3, NUM_NORMAL_STATES, 3, (Q+1)+(T+1) );
      MATRIX_3D_Fill( test_MX3, -INF );
      EDGEBOUNDS_Reuse( test_edg, Q, T );
      MATRIX_2D_Fill( sp_MX, -INF );

      num_writes = 0;
      num_clears = 0;
   }
   #endif

   /* --------------------------------------------------------------------------------- */

   /* initialize logsum lookup table if it has not already been */
   logsum_Init();

   /* query sequence */
   mx          = st_SMX;
   seq         = query->seq;
   N           = edg->N;
   /* local or global alignments? */
   is_local    = target->isLocal;
   sc_E        = (is_local) ? 0 : -INF;

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

   /* initialize indexes */
   q_0   = 0;
   r_0b  = 0;
   r_0e  = 0;

   /* UNROLLED INITIAL ROW */
   {
      /* check if query position is in domain */
      is_q_0_in_dom_range = IS_IN_RANGE( Q_range.beg, Q_range.end, q_0 );
      /* get edgebound range */
      EDGEBOUNDS_NxtRow( edg, &r_0b, &r_0e, q_0 );

      /* initialize special states */
      XMX(SP_E, q_0) = -INF; 
      XMX(SP_J, q_0) = -INF; 
      XMX(SP_C, q_0) = -INF; 
      /* S->N, p=1 */
      XMX(SP_N, q_0) = 0.0f; 
      /* S->N->B, no N-tail */
      XMX(SP_B, q_0) = XSC(SP_N,SP_MOVE); 

      /* only compute if in domain range */
      if ( is_q_0_in_dom_range == true )
      {
         /* FOR every BOUND in zero row (-INF values set during initialization, so unneccessary) */
         for (r_0 = r_0b; r_0 < r_0e; r_0++)
         {
            /* get bound data */
            bnd   = &EDG_X(edg, r_0);
            id    = bnd->id;
            lb_T  = bnd->lb <= 0;
            lb_0  = MAX(bnd->lb, T_range.beg);    /* can't overflow the left edge */
            rb_T  = bnd->rb >= T;
            rb_0  = MIN(bnd->rb, T_range.end);    /* can't overflow the right edge */

            /* fetch data mapping bound start location to data block in sparse matrix */
            qx0 = VECTOR_INT_Get( st_SMX->imap_cur, r_0 );    /* (q_0, t_0) location offset */
            qx1 = VECTOR_INT_Get( st_SMX->imap_prv, r_0 );    /* (q_1, t_0) location offset */

            /* initial location for square matrix and mapping to sparse matrix */
            t_0 = lb_0;
            tx0 = t_0 - bnd->lb;    /* total_offset = offset_location - starting_location */

            /* FOR every position in TARGET profile */
            for (t_0 = lb_0; t_0 < rb_0; t_0++)
            {
               tx0 = t_0 - bnd->lb;
               /* zero column is -inf in logspace.  We can skip this step and convert to normal space now. */
               MSMX(qx0, tx0) = -INF;
               ISMX(qx0, tx0) = -INF;
               DSMX(qx0, tx0) = -INF; 
            }
         }
      }

      /* lookback one row */
      r_1b = r_0b;
      r_1e = r_0e;
   }
   
   /* MAIN RECURSION */
   /* FOR every position in QUERY sequence (row in matrix) */
   for (q_0 = 1; q_0 <= Q; q_0++)
   {
      q_1 = q_0 - 1;
      t_0 = 0;

      /* check if query position is in domain */
      is_q_0_in_dom_range = (q_0 >= Q_range.beg && q_0 < Q_range.end);
      /* get edgebound range */
      EDGEBOUNDS_NxtRow( edg, &r_0b, &r_0e, q_0 );

      /* Get next sequence character */
      a = seq[q_1];  /* off-by-one */
      A = AA_REV[a];

      XMX(SP_E, q_0) = -INF;

      /* only compute if in domain range */
      if ( is_q_0_in_dom_range == true )
      {
         /* FOR every BOUND in current ROW */
         for (r_0 = r_0b; r_0 < r_0e; r_0++)
         {
            /* get bound data */
            bnd   = &EDG_X(edg, r_0);
            id    = bnd->id;
            lb_T  = bnd->lb <= 0;
            lb_0  = MAX(bnd->lb, T_range.beg);   /* can't overflow left edge */
            rb_T  = bnd->rb >= T;
            rb_0  = MIN(bnd->rb, T_range.end);   /* can't overflow right edge */

            /* fetch data mapping bound start location to data block in sparse matrix */
            qx0 = VECTOR_INT_Get( st_SMX->imap_cur, r_0 );    /* (q_0, t_0) location offset */
            qx1 = VECTOR_INT_Get( st_SMX->imap_prv, r_0 );    /* (q_1, t_0) location offset */

            /* initial location for square matrix and mapping to sparse matrix */
            t_0 = lb_0;
            tx0 = t_0 - bnd->lb;    /* total_offset = offset_location - starting_location */
            tx1 = tx0 - 1;

            /* unrolled first loop: special case for left edge of range */
            if ( true )
            {
               t_0 = lb_0;
               tx0 = t_0 - bnd->lb;

               /* zero column is -inf in logspace.  We can skip this step and convert to normal space now. */
               MSMX(qx0, tx0) = -INF;
               ISMX(qx0, tx0) = -INF;
               DSMX(qx0, tx0) = -INF; 
            }

            /* MAIN RECURSION */
            /* FOR every position in TARGET profile */
            for (t_0 = lb_0 + 1; t_0 < rb_0 - 1; t_0++)
            {
               t_1 = t_0 - 1; 
               tx0 = t_0 - bnd->lb;
               tx1 = tx0 - 1;

               /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
               /* best previous state transition (match takes the diag element of each prev state) */
               prv_M = MSMX(qx1, tx1) + TSC(t_1, M2M);
               prv_I = ISMX(qx1, tx1) + TSC(t_1, I2M);
               prv_D = DSMX(qx1, tx1) + TSC(t_1, D2M);
               prv_B = XMX(SP_B, q_1) + TSC(t_1, B2M); /* from begin match state (new alignment) */
               /* best-to-match */
               prv_sum = logsum( logsum( prv_M, prv_I ),
                                 logsum( prv_B, prv_D ) );
               MSMX(qx0, tx0) = prv_sum + MSC(t_0, A);

               /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
               /* previous states (match takes the previous row (upper) of each state) */
               prv_M = MSMX(qx1, tx0) + TSC(t_0, M2I);
               prv_I = ISMX(qx1, tx0) + TSC(t_0, I2I);
               /* best-to-insert */
               prv_sum = logsum( prv_M, prv_I );
               ISMX(qx0, tx0) = prv_sum + ISC(t_0, A);

               // if (q_0 == 1 && t_0 > 12 && t_0 < 15 ) {
               //    printf("(q_0,t_0)=(%d,%d): MMX=%f, IMX=%f => IMX=%f\n", 
               //       q_0, t_0, MSMX(qx1, tx0), ISMX(qx1, tx0), ISMX(qx0, tx0));
               // }

               /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) */
               /* previous states (match takes the previous column (left) of each state) */
               prv_M = MSMX(qx0, tx1) + TSC(t_1, M2D);
               prv_D = DSMX(qx0, tx1) + TSC(t_1, D2D);
               /* best-to-delete */
               prv_sum = logsum( prv_M, prv_D );
               DSMX(qx0, tx0) = prv_sum;


               /* UPDATE E STATE */
               prv_M = MSMX(qx0, tx0) + sc_E;
               prv_D = DSMX(qx0, tx0) + sc_E;
               prv_E = XMX(SP_E, q_0);
               /* best-to-e-state */
               prv_sum = logsum( logsum(  prv_M, prv_D ),
                                          prv_E );
               XMX(SP_E, q_0) = prv_sum;

               /* embed linear row into quadratic test matrix */
               #if DEBUG
               {
                  MX_2D(cloud_MX, q_0, t_0) = 1.0;
                  MX_3D(test_MX, MAT_ST, q_0, t_0) = MSMX(qx0, tx0);
                  MX_3D(test_MX, INS_ST, q_0, t_0) = ISMX(qx0, tx0);
                  MX_3D(test_MX, DEL_ST, q_0, t_0) = DSMX(qx0, tx0);
               }
               #endif
            }

            /* unrolled final loop: special case for right edge of range */
            if (true && rb_0 > 1)  
            {
               // printf("q_0,t_0 = %d,%d\n", q_0, t_0);
               t_0 = rb_0 - 1;
               t_1 = t_0 - 1; 
               tx0 = t_0 - bnd->lb;
               tx1 = tx0 - 1;

               /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
               /* best previous state transition (match takes the diag element of each prev state) */
               prv_M = MSMX(qx1, tx1) + TSC(t_1, M2M);
               prv_I = ISMX(qx1, tx1) + TSC(t_1, I2M);
               prv_D = DSMX(qx1, tx1) + TSC(t_1, D2M);
               prv_B = XMX(SP_B, q_1) + TSC(t_1, B2M);    /* from begin match state (new alignment) */
               /* sum-to-match */
               prv_sum = logsum( logsum( prv_M, prv_I ),
                                 logsum( prv_D, prv_B ) );
               MSMX(qx0, tx0) = prv_sum + MSC(t_0, A);

               /* FIND SUM OF PATHS TO INSERT STATE (unrolled) */
               ISMX(qx0, tx0) = -INF;

               /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) (unrolled) */
               /* previous states (match takes the left element of each state) */
               prv_M = MSMX(qx0, tx1) + TSC(t_1, M2D);
               prv_D = DSMX(qx0, tx1) + TSC(t_1, D2D);
               /* sum-to-delete */
               prv_sum = logsum( prv_M, prv_D );
               DSMX(qx0, tx0) = prv_sum;

               /* UPDATE E STATE (unrolled) */
               prv_E = XMX(SP_E, q_0);
               prv_M = MSMX(qx0, tx0);
               prv_D = DSMX(qx0, tx0);
               /* best-to-begin */
               prv_sum = logsum( logsum(  prv_D, prv_M ),
                                          prv_E ); 
               XMX(SP_E, q_0) = prv_sum;

               /* embed linear row into quadratic test matrix */
               #if DEBUG
               {
                  MX_2D(cloud_MX, q_0, t_0) = 1.0;
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
      prv_J = XMX(SP_J, q_1) + XSC(SP_J, SP_LOOP);       /* J->J */
      prv_E = XMX(SP_E, q_0) + XSC(SP_E, SP_LOOP);       /* E->J is E's "loop" */
      prv_sum = logsum( prv_J, prv_E );
      XMX(SP_J, q_0) = prv_sum;

      /* C state */
      prv_C = XMX(SP_C, q_1) + XSC(SP_C, SP_LOOP);
      prv_E = XMX(SP_E, q_0) + XSC(SP_E, SP_MOVE);
      prv_sum = logsum( prv_C, prv_E );
      XMX(SP_C, q_0) = prv_sum;

      /* N state */
      prv_N = XMX(SP_N, q_1) + XSC(SP_N, SP_LOOP);
      XMX(SP_N, q_0) = prv_N;

      /* B state */
      prv_N = XMX(SP_N, q_0) + XSC(SP_N, SP_MOVE);         /* N->B is N's move */
      prv_J = XMX(SP_J, q_0) + XSC(SP_J, SP_MOVE);         /* J->B is J's move */
      prv_sum = logsum( prv_N, prv_J ); 
      XMX(SP_B, q_0) = prv_sum;

      /* SET CURRENT ROW TO PREVIOUS ROW */
      r_1b = r_0b;
      r_1e = r_0e;
   }

   /* T state */
   sc_best     = XMX(SP_C, Q) + XSC(SP_C, SP_MOVE);
   *sc_final   = sc_best;

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
int run_Bound_Backward_Sparse (  const SEQUENCE*      query,         /* query sequence */
                                 const HMM_PROFILE*   target,        /* target HMM model */
                                 const int            Q,             /* query length */
                                 const int            T,             /* target length */
                                 MATRIX_3D_SPARSE*    st_SMX,        /* normal state matrix */
                                 MATRIX_2D*           sp_MX,         /* special state matrix */
                                 EDGEBOUNDS*          edg,           /* edgebounds */
                                 RANGE*               dom_range,     /* (OPTIONAL) domain range for computing fwd/bck on specific domain. If NULL, computes complete fwd/bck. */
                                 float*               sc_final )     /* (OUTPUT) final score */
{
   /* vars for matrix access */
   MATRIX_3D_SPARSE*    n_mx;                /* normal state matrix */
   MATRIX_2D*           s_mx;                /* special state matrix */

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
   BOUND*   bnd;                             /* current bound */
   int      id;                              /* id in edgebound list (row/diag) */
   int      r_0;                             /* current index for current row */
   int      r_0b, r_0e;                      /* begin and end indices for current row in edgebound list */
   int      r_1;                             /* current index for previous row */
   int      r_1b, r_1e;                      /* begin and end indices for current row in edgebound list */
   int      le_0, re_0;                      /* right/left matrix bounds of current diag */
   int      lb_0, rb_0;                      /* bounds of current search space on current diag */
   int      lb_1, rb_1;                      /* bounds of current search space on previous diag */
   int      lb_2, rb_2;                      /* bounds of current search space on 2-back diag */
   bool     lb_T, rb_T;                      /* checks if edge touches right bound of matrix */

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
   logsum_Init();

   /* query sequence */
   mx          = st_SMX;
   seq         = query->seq;
   N           = edg->N;
   /* local or global alignments? */
   is_local    = target->isLocal;
   sc_E        = (is_local) ? 0 : -INF;

   /* domain range (query sequence) */
   if (dom_range == NULL) {
      Q_range.beg = 0;
      Q_range.end = Q;
   } else {
      Q_range = *dom_range;
   }
   /* target range */
   T_range.beg = 1;
   T_range.end = T;

   /* init indexes */
   q_0   = Q;
   r_0b  = edg->N - 1;
   r_0e  = edg->N - 1;

   /* UNROLLED INITIAL Qth ROW */
   {
      /* if inside domain */
      is_q_0_in_dom_range = IS_IN_RANGE( Q_range.beg, Q_range.end, q_0 );
      /* get edgebound range */
      EDGEBOUNDS_PrvRow( edg, &r_0b, &r_0e, q_0 );

      /* intialize special states */
      XMX(SP_J, q_0) = -INF;
      XMX(SP_B, q_0) = -INF;
      XMX(SP_N, q_0) = -INF;
      XMX(SP_C, q_0) = XSC(SP_C, SP_MOVE);
      XMX(SP_E, q_0) = XMX(SP_C, q_0) + XSC(SP_E, SP_MOVE);

      /* if sequence position is in domain range */ 
      if ( is_q_0_in_dom_range == true )
      {
         /* if Q-row bounds are not empty and the right-most bound spans T ( covers bottom-right corner ) */
         if ( (r_0b - r_0e > 0) && (EDG_X(edg, r_0b).rb > T) )
         {
            /* get bound data */
            bnd   = &EDG_X(edg, r_0b);
            /* fetch data location to bound start location (in offset) */
            qx0   = VECTOR_INT_Get( st_SMX->imap_cur, r_0b );     /* (q_0, t_0) location offset */

            /* location for square matrix and mapping to sparse matrix */
            t_0 = T;
            tx0 = t_0 - bnd->lb;    /* total_offset = offset_location - starting_location */

            MSMX(qx0, tx0) = XMX(SP_E, q_0);
            ISMX(qx0, tx0) = -INF;
            DSMX(qx0, tx0) = XMX(SP_E, q_0);

            #if DEBUG 
            {
               MX_2D(cloud_MX, q_0, t_0) = 1.0;
               MX_3D(test_MX, MAT_ST, q_0, t_0) = MSMX(qx0, tx0);
               MX_3D(test_MX, INS_ST, q_0, t_0) = ISMX(qx0, tx0);
               MX_3D(test_MX, DEL_ST, q_0, t_0) = DSMX(qx0, tx0);
            }
            #endif
         }

         /* Initialize normal states (sparse) */
         for (r_0 = r_0b; r_0 > r_0e; r_0--) 
         {
            /* get bound data */
            bnd   = &EDG_X(edg, r_0);            /* bounds for current bound */
            id    = bnd->id;
            lb_T  = bnd->lb <= 0;
            lb_0  = MAX(bnd->lb, T_range.beg);   /* can't overflow left edge */
            rb_T  = bnd->rb >= T;
            rb_0  = MIN(bnd->rb, T_range.end);   /* can't overflow right edge */
            printf("INIT q=%d, r=%d, {%d,%d,%d}\n", q_0, r_0, id, bnd->lb, bnd->rb);

            /* fetch data mapping bound start location to data block in sparse matrix */
            qx0   = VECTOR_INT_Get( st_SMX->imap_cur, r_0 );    /* (q_0, t_0) location offset */
            qx1   = VECTOR_INT_Get( st_SMX->imap_nxt, r_0 );    /* (q_1, t_0) location offset */

            /* location for square matrix and mapping to sparse matrix */
            t_0 = lb_0;
            tx0 = t_0 - bnd->lb;    /* total_offset = offset_location - starting_location */

            for (t_0 = rb_0 - 1; t_0 >= lb_0; t_0--)
            {
               /* real target index */
               t_1 = t_0 + 1;
               /* calculate offset from beginning of sparse data block */ 
               tx0 = t_0 - bnd->lb;
               tx1 = tx0 + 1;

               prv_E = XMX(SP_E, Q) + sc_E;
               prv_D = DSMX(qx0, tx1)  + TSC(t_0, M2D);
               MSMX(qx0, tx0) = logsum( prv_E, prv_D );

               prv_E = XMX(SP_E, Q) + sc_E;
               prv_D = DSMX(qx0, tx1)  + TSC(t_0, D2D);
               DSMX(qx0, tx0) = logsum( prv_E, prv_D );

               ISMX(qx0, tx0) = -INF;

               #if DEBUG 
               {
                  MX_2D(cloud_MX, q_0, t_0) = 1.0;
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
   

   /* MAIN RECURSION */
   /* FOR every bound in EDGEBOUND */
   for (q_0 = Q-1; q_0 > 0; q_0--)
   {
      q_1 = q_0 + 1;
      t_0 = 0;

      /* if inside domain */
      is_q_0_in_dom_range = IS_IN_RANGE( Q_range.beg, Q_range.end, q_0 );
      is_q_1_in_dom_range = IS_IN_RANGE( Q_range.beg, Q_range.end, q_1 );
      /* get edgebound range */
      EDGEBOUNDS_PrvRow( edg, &r_0b, &r_0e, q_0 );

      /* Get next sequence character */
      a = seq[q_0];
      A = AA_REV[a];

      /* init B STATE (sparse) */
      XMX(SP_B, q_0) = -INF;
      /* if previous q is in domain range, update B state */
      if ( is_q_1_in_dom_range == true )
      {
         for (r_1 = r_1b; r_1 > r_1e; r_1--) 
         {
            /* get bound data */
            bnd   = &EDG_X(edg, r_1);              /* bounds for current bound */
            id    = bnd->id;
            lb_T  = bnd->lb <= 0;
            lb_0  = MAX(bnd->lb, T_range.beg);   /* can't overflow left edge */
            rb_T  = bnd->rb >= T;
            rb_0  = MIN(bnd->rb, T_range.end);   /* can't overflow right edge */

            /* fetch data location to bound start location (in offset) */
            qx1 = VECTOR_INT_Get( st_SMX->imap_cur, r_1 );    /* (q_0, t_0) location offset */

            /* location for square matrix and mapping to sparse matrix */
            t_0 = lb_0;
            tx0 = t_0 - bnd->lb;    /* total_offset = offset_location - starting_location */

            for (t_0 = rb_0 - 1; t_0 >= lb_0; t_0--)
            {
               t_1 = t_0 - 1;
               /* calculate offset from beginning of sparse data block */ 
               tx0 = t_0 - bnd->lb;
               tx1 = tx0 - 1;

               prv_sum = XMX(SP_B, q_0);
               prv_M = MSMX(qx1, tx0) + TSC(t_1, B2M) + MSC(t_0, A);
               XMX(SP_B, q_0) = logsum( prv_sum, prv_M);
            }
         }
      }

      prv_J = XMX(SP_J, q_1) + XSC(SP_J, SP_LOOP);
      prv_B = XMX(SP_B, q_0) + XSC(SP_J, SP_MOVE);
      XMX(SP_J, q_0) = logsum( prv_J, prv_B );

      prv_C = XMX(SP_C, q_1) + XSC(SP_C, SP_LOOP);
      XMX(SP_C, q_0) = prv_C;

      prv_J = XMX(SP_J, q_0) + XSC(SP_E, SP_LOOP);
      prv_C = XMX(SP_C, q_0) + XSC(SP_E, SP_MOVE);
      XMX(SP_E, q_0) = logsum( prv_J, prv_C );

      prv_N = XMX(SP_N, q_1) + XSC(SP_N, SP_LOOP);
      prv_B = XMX(SP_B, q_0) + XSC(SP_N, SP_MOVE);
      XMX(SP_N, q_0) = logsum( prv_N, prv_B );

      if ( is_q_0_in_dom_range == true )
      {
         /* if there is a bound on row and the right-most bound spans T (right edge of matrix) */
         if ( (r_0b - r_0e > 0) && (EDG_X(edg, r_0b).rb > T) )
         {
            /* get bound data */
            bnd   = &EDG_X(edg, r_0b);       /* bounds for current bound */
            id    = bnd->id;
            lb_T  = bnd->lb <= 0;
            lb_0  = MAX(bnd->lb, T_range.beg);   /* can't overflow left edge */
            rb_T  = bnd->rb >= T;
            rb_0  = MIN(bnd->rb, T_range.end);   /* can't overflow right edge */

            /* fetch data location to bound start location (in offset) */
            qx0 = VECTOR_INT_Get( st_SMX->imap_cur, r_0b );    /* (q_0, t_0) location offset */

            /* location for square matrix and mapping to sparse matrix */
            t_0 = T;
            tx0 = t_0 - bnd->lb;    /* total_offset = offset_location - starting_location */

            MSMX(qx0, tx0) = XMX(SP_E, q_0);
            ISMX(qx0, tx0) = -INF;
            DSMX(qx0, tx0) = XMX(SP_E, q_0);

            #if DEBUG 
            {
               MX_2D(cloud_MX, q_0, t_0) = 1.0;
               MX_3D(test_MX, MAT_ST, q_0, t_0) = MSMX(qx0, tx0);
               MX_3D(test_MX, INS_ST, q_0, t_0) = ISMX(qx0, tx0);
               MX_3D(test_MX, DEL_ST, q_0, t_0) = DSMX(qx0, tx0);
            }
            #endif
         }

         /* FOR every EDGEBOUND in current ROW */
         for (r_0 = r_0b; r_0 > r_0e; r_0--)
         {
            /* get bound data */
            bnd   = &EDG_X(edg, r_0);              /* bounds for current bound */
            id    = bnd->id;
            lb_T  = bnd->lb <= 0;
            lb_0  = MAX(bnd->lb, T_range.beg);   /* can't overflow left edge */
            rb_T  = bnd->rb >= T;
            rb_0  = MIN(bnd->rb, T_range.end);   /* can't overflow right edge */
            
            /* fetch data location to bound start location (in offset) */
            /* NOTE */
            qx0   = VECTOR_INT_Get( st_SMX->imap_cur, r_0 );    /* (q_0, t_0) location offset */
            qx1   = VECTOR_INT_Get( st_SMX->imap_nxt, r_0 );    /* (q_0, t_0) location offset */

            /* location for square matrix and mapping to sparse matrix */
            t_0 = T;
            tx0 = t_0 - bnd->lb;    /* total_offset = offset_location - starting_location */

            /* FOR every position in TARGET profile */
            for (t_0 = rb_0 - 1; t_0 >= lb_0; t_0--)
            {
               t_1 = t_0 + 1;
               /* calculate offset from beginning of sparse data block */ 
               tx0 = t_0 - bnd->lb;
               tx1 = tx0 + 1;

               /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
               prv_M = MSMX(qx1, tx1) + TSC(t_0, M2M) + MSC(t_1, A);
               prv_I = ISMX(qx1, tx0) + TSC(t_0, M2I) + ISC(t_1, A);
               prv_D = DSMX(qx0, tx1) + TSC(t_0, M2D);
               prv_E = XMX(SP_E, q_0) + sc_E;     /* from end match state (new alignment) */
               /* best-to-match */
               prv_sum = logsum( 
                              logsum( prv_M, prv_I ),
                              logsum( prv_E, prv_D ) );
               MSMX(qx0, tx0) = prv_sum;

               /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
               prv_M = MSMX(qx1, tx1) + TSC(t_0, I2M) + MSC(t_1, A);
               prv_I = ISMX(qx1, tx0) + TSC(t_0, I2I) + ISC(t_0, A);
               /* best-to-insert */
               prv_sum = logsum( prv_M, prv_I );
               ISMX(qx0, tx0) = prv_sum;

               /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
               prv_M = MSMX(qx1, tx1) + TSC(t_0, D2M) + MSC(t_1, A);
               prv_D = DSMX(qx0, tx1) + TSC(t_0, D2D);
               prv_E = XMX(SP_E, q_0) + sc_E;
               /* best-to-delete */
               prv_sum = logsum( prv_M, 
                              logsum( prv_D, prv_E ) );
               DSMX(qx0, tx0) = prv_sum;

               #if DEBUG 
               {
                  MX_2D(cloud_MX, q_0, t_0) = 1.0;
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

   /* FINAL ROW (i = 0) */
   /* At q_0 = 0, only N,B states are reachable. */
   q_0 = 0;
   q_1 = q_0 + 1;

   t_0 = 0;
   t_1 = t_0 + 1;

   a = seq[q_0];
   A = AA_REV[a];

   /* get edgebound range */
   EDGEBOUNDS_PrvRow( edg, &r_0b, &r_0e, q_0 );

   /* FINAL i = 0 row */
   a = seq[0];
   A = AA_REV[a];

   /* B STATE (SPARSE) */
   XMX(SP_B, q_0) = -INF;
   /* if previous q is in domain, update B state */
   if ( is_q_1_in_dom_range == true )
   {
      for (r_1 = r_1b; r_1 > r_1e; r_1--) 
      {
         bnd   = &EDG_X(edg, r_1);         /* bounds for current bound */
         id    = bnd->id;
         lb_T  = bnd->lb <= 0;
         lb_0  = MAX(bnd->lb, T_range.beg);   /* can't overflow left edge */
         rb_T  = bnd->rb >= T;
         rb_0  = MIN(bnd->rb, T_range.end);   /* can't overflow right edge */

         /* fetch data location to bound start location (in offset) */
         qx1 = VECTOR_INT_Get( st_SMX->imap_cur, r_1 );    /* (q_0, t_0) location offset */

         /* location for square matrix and mapping to sparse matrix */
         t_0 = lb_0;
         tx0 = t_0 - bnd->lb;    /* total_offset = offset_location - starting_location */

         for (t_0 = rb_0 - 1; t_0 >= lb_0; t_0--)
         {
            t_1 = t_0 - 1;
            /* calculate offset from beginning of sparse data block */ 
            tx0 = t_0 - bnd->lb;
            tx1 = tx0 - 1;

            prv_sum = XMX(SP_B, q_0);
            prv_M = MSMX(qx1, tx0) + TSC(t_1, B2M) + MSC(t_0, A);
            XMX(SP_B, q_0) = logsum( prv_sum, prv_M );
         }
      }
   }

   XMX(SP_J, q_0) = -INF;
   XMX(SP_C, q_0) = -INF;
   XMX(SP_E, q_0) = -INF;

   prv_N = XMX(SP_N, q_1) + XSC(SP_N, SP_LOOP);
   prv_B = XMX(SP_B, q_0) + XSC(SP_N, SP_MOVE);
   XMX(SP_N, q_0) = logsum( prv_N, prv_B );

   sc_best     = XMX(SP_N, 0);
   *sc_final   = sc_best;

   return STATUS_SUCCESS;
}
