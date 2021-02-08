/*******************************************************************************
 *  FILE:      posterior_traceback_sparse.h
 *  PURPOSE:   The Maximum Posterior Probability and Optimal Alignment.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
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

/* header */
#include "_algs_sparse.h"
#include "posterior_traceback_sparse.h"

/* This determines whether optimal alignment will be computed for cell, prevent NaN errors */
#define TSC_DELTA(t_0,tr)     ((TSC((t_0),(tr)) == -INF) ? FLT_MIN : 1.0)
#define XSC_DELTA(sp,tr)      ((XSC((sp),(tr)) == -INF) ? FLT_MIN : 1.0)

/*! FUNCTION:  run_Posterior_Traceback_Sparse()
 *  SYNOPSIS:  Traceback of Posterior Probability matrix <...post>.
 *             As there is no transition probabilities, this simply backtraces 
 *             by taking the greedy algorithm and choosing the maximal scoring next state
 *             in the graph.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_Posterior_Optimal_Traceback_Sparse(     const SEQUENCE*         query,            /* query sequence */
                                            const HMM_PROFILE*      target,           /* target hmm model */
                                            const int               Q,                /* query length */
                                            const int               T,                /* target length */
                                            EDGEBOUNDS*             edg,              /* edgebounds */
                                            RANGE*                  dom_range,        /* query span of bounds */
                                            MATRIX_3D_SPARSE*       st_SMX_opt,       /* optimal accuracy normal matrix */
                                            MATRIX_2D*              sp_MX_opt,        /* optimal accuracy special matrix */
                                            ALIGNMENT*              aln )             /* OUTPUT: optimal alignment */
{
   printf("=== run_Posterior_Traceback_Sparse() [BEGIN] ===\n");
   FILE*                fp;

   /* generic dp matrix pointers */
   MATRIX_3D_SPARSE*    st_SMX;
   MATRIX_2D*           sp_MX;

   /* vars for accessing query/target data structs */
   char     a;                               /* store current character in sequence */
   int      A;                               /* store int value of character */
   char*    seq;                             /* alias for getting seq */
   bool     is_local;                        /* whether local or global alignment */

   /* vars for indexing into data matrices by row-col */
   int      b, d, i, j, k;                   /* antidiagonal, row, column indices */
   int      q_0, q_1;                        /* real index of current and previous rows (query) */
   int      qx0, qx1;                        /* maps column index into data index (query) */
   int      t_0, t_1;                        /* real index of current and previous columns (target) */
   int      tx0, tx1;                        /* maps target index into data index (target)  */
   int      t_range;                         /* range of targets on current row */

   /* vars for indexing into edgebound lists */
   BOUND*   bnd;                             /* current bound */
   int      id_0;                            /* id in edgebound list (row/diag) */
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
   float    cur;                             /* current state */
   float    prv_M, prv_I, prv_D;             /* previous (M) match, (I) insert, (D) delete states */
   float    prv_B, prv_E;                    /* previous (B) begin and (E) end states */
   float    prv_J, prv_N, prv_C;             /* previous (J) jump, (N) initial, and (C) terminal states */
   float    prv_loop, prv_move;              /* previous loop and move for special states */
   float    prv_sum, prv_best;               /* temp subtotaling vars */
   float    sc_best;                         /* final best scores */
   float    sc_M, sc_I, sc_D, sc_E;          /* match, insert, delete, end scores */

   /* vars for alignment traceback */
   float          sc_cur, sc_max;            /* current and maximum score */
   int            st_cur, st_prv, st_max;    /* current, previous state in traceback */
   int            q_prv, t_prv, t_max;       /* previous and maximum query and target positions */
   TRACE*         tr;                        /* trace object for appending */
   const float    tol = 1e-5;                /* acceptable tolerance range for "equality tests" */
   bool           is_found;                  /* determines whether cell is found */
   int            aln_len;                   /* length of alignment */

   /* --------------------------------------------------------------------------- */

   /* open test file */
   #if DEBUG
   {
      fp = fopen("test_output/my.optimal_traceback.tsv", "w+");
   }
   #endif

   /* get sequence and trace */
   st_SMX   = st_SMX_opt;
   sp_MX    = sp_MX_opt;
   seq      = query->seq; 
   tr       = aln->traces->data; 
   /* local or global? */
   is_local = target->isLocal;

   /* intial indexes */
   q_0   = q_prv = Q;
   t_0   = t_prv = 0;
   r_0b  = edg->N - 1;
   r_0e  = edg->N - 1;
   /* get edgebound range */
   EDGEBOUNDS_PrvRow( edg, &r_0b, &r_0e, q_0 );

   /* clear memory for trace */
   ALIGNMENT_Reuse( aln, Q, T );

   /* Backtracing, so C is end state */
   ALIGNMENT_Append( aln, T_ST, q_0, t_0 );
   ALIGNMENT_Append( aln, C_ST, q_0, t_0 );
   st_prv = C_ST;

   /* End of trace is S state */
   while (st_prv != S_ST)
   {
      /* if we have reached the end of the query sequence */
      if (q_0 == 0) {
         ALIGNMENT_Append( aln, S_ST, q_0, t_0 );
         break;
      } 

      /* if q_0 has been decremented, then get edgebound range of next row  */
      if ( q_prv != q_0 ) {
         /* get edgebound range */
         EDGEBOUNDS_PrvRow( edg, &r_0b, &r_0e, q_0 );
      }

      /*! TODO: (?) can be optimized by only updating during certain cases */
      is_found = false;
      /* if in the core model, find location */
      if ( st_prv == M_ST || st_prv == I_ST || st_prv == D_ST )
      {
         /* if either t_0 or q_0 have been decremented, then lookup index offset */
         if ( t_prv != t_0 || q_prv != q_0 ) {
            /* search through bounds in row for bound containing index */
            for ( r_0 = r_0b; r_0 > r_0e; r_0-- ) {
               /* find bound that contains current index */
               bnd = &EDG_X(edg, r_0);
               /* if t_0 is inside range of current bound */
               if ( t_0 >= bnd->lb && t_0 < bnd->rb ) {
                  /* fetch data mapping bound start location to data block in sparse matrix */
                  qx0 = VECTOR_INT_Get( st_SMX->imap_cur, r_0 );    /* (q_0, 0) location offset */
                  qx1 = VECTOR_INT_Get( st_SMX->imap_prv, r_0 );    /* (q_1, 0) location offset */
                  /* total_offset = offset_location - starting_location */
                  tx0 = t_0 - bnd->lb;    
                  tx1 = tx0 - 1;
                  /* found location, so break from loop */
                  is_found = true;
                  break;
               }
            }
         }
         /* if not found, throw error */
         if ( is_found == false ) {
            fprintf( stderr, "ERROR: Impossible position in model reached at {%s,%d,%d}\n", 
               STATE_NAMES[st_prv], q_0, t_0);
            exit(EXIT_FAILURE);
         }
      } 

      /* for finding next maximal state, init max to -inf */
      sc_max = -INF;

      /* update previous */
      q_prv = q_0;
      t_prv = t_0; 

      /* previous target and query sequence */
      q_1 = q_0 - 1;
      t_1 = t_0 - 1;
      
      /* get next sequence character */
      a = seq[q_1];
      A = AA_REV[a];

      /* jump from current state to the prev state */
      switch ( st_prv )
      {
         /* C STATE to {C,E} */
         case C_ST:  /* C(q_0) comes from C(q_1) or E(q_0) */
         {
            /* current state */
            sc_cur = XMX(SP_C, q_0);
            /* No valid alignment goes to -INF */
            if ( sc_cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible C_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }
            
            /* possible previous states */
            prv_C = XMX(SP_C, q_1);
            prv_E = XMX(SP_E, q_0);

            if ( sc_max < prv_C ) {
               st_cur = C_ST;
               sc_max = prv_C;
            }
            if ( sc_max < prv_E ) {
               st_cur = E_ST;
               sc_max = prv_E;
            }
            if ( sc_max == -INF ) {
               fprintf( stderr, "ERROR: Failed to trace from B_ST at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }
         } break;

         /* E STATE to {M, D} */
         case E_ST:  /* E connects from any M state. t_0 is set here. */
         {
            /* current state */
            sc_cur = XMX(SP_E, q_0);
            /* check for error state */
            if ( sc_cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible E_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }

            if ( is_local )  /* local mode: ends in M */
            {
               /* can't come from D, in a *local* Viterbi alignment. */
               st_cur = M_ST;

               /* possible previous states (any M state) */
               /* FOR every BOUND in current ROW */
               for ( r_0 = r_0b; r_0 > r_0e; r_0-- )
               {
                  /* get bound data */
                  bnd   = &EDG_X(edg, r_0);
                  id_0  = bnd->id;
                  lb_0  = bnd->lb;
                  rb_0  = bnd->rb;

                  /* fetch data mapping bound start location to data block in sparse matrix */
                  qx0 = VECTOR_INT_Get( st_SMX->imap_cur, r_0 );    /* (q_0, t_0) location offset */
                  qx1 = VECTOR_INT_Get( st_SMX->imap_prv, r_0 );    /* (q_1, t_0) location offset */

                  /* initial location for square matrix and mapping to sparse matrix */
                  t_0 = lb_0;
                  tx0 = t_0 - bnd->lb;    /* total_offset = offset_location - starting_location */

                  for ( t_0 = lb_0; t_0 < rb_0; t_0++ ) 
                  {
                     t_1 = t_0 - 1; 
                     tx0 = t_0 - bnd->lb;
                     tx1 = tx0 - 1;

                     /* possible previous state */
                     prv_M = MSMX(qx0, tx0);

                     /* check if the maximal scoring state found */
                     if ( sc_max < prv_M ) {
                        sc_max = prv_M;
                        t_max = t_0;
                     }
                  }
               }
               t_0 = t_max;
            }
            else /*! TODO: (?) glocal mode: we either come from D_M or M_M */
            {
               fprintf(stderr, "ERROR: Glocal mode not supported for sparse alignment.");
               ERRORCHECK_exit(STATUS_FAILURE);
            }
         } break;

         /* M STATE to {B,M,I,D} */
         case M_ST:  /* M connects from (q_1, t_1), or B */
         {
            /* current score */
            sc_cur = MSMX(qx0, tx0);
            /* No valid alignment goes to -INF */
            if ( sc_cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible M_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }

            /* possible previous states */
            prv_B = XMX(SP_B, q_1);
            prv_M = MSMX(qx1, tx1);
            prv_I = ISMX(qx1, tx1);
            prv_D = DSMX(qx1, tx1);

            /* find maximum next state score */
            if ( sc_max < prv_B ) {
               st_cur = B_ST;
               sc_max = prv_B;
            }
            if ( sc_max < prv_M ) {
               st_cur = M_ST;
               sc_max = prv_M;
            }
            if ( sc_max < prv_I ) {
               st_cur = I_ST;
               sc_max = prv_I;
            }
            if ( sc_max < prv_D ) {
               st_cur = D_ST;
               sc_max = prv_D;
            }
            if ( sc_max == -INF ) {
               fprintf( stderr, "ERROR: Failed to trace from M_ST at (%d,%d)\n", t_0, q_0);
               exit(EXIT_FAILURE);
            }

            /* update index to previous state */
            t_0--; q_0--;
         } break;

         /* D STATE to {M,D} */
         case D_ST:  /* D connects from M,D at (q_0, t_1) */
         {
            /* current state */
            sc_cur = DSMX(qx0, tx0);
            t_1 = t_0 - 1;
            /* No valid alignment goes to -INF */
            if ( sc_cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible D_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }

            /* possible previous states */
            prv_M = MSMX(qx0, tx1);
            prv_D = DSMX(qx0, tx1);

            /* find maximum next state score */            
            if ( sc_max < prv_M ) {
               st_cur = M_ST;
               sc_max = prv_M;
            }
            if ( sc_max < prv_D ) {
               st_cur = D_ST;
               sc_max = prv_D;
            }
            if ( sc_max == -INF ) {
               fprintf( stderr, "ERROR: Failed to trace from D_ST at (%d,%d)\n", t_0, q_0);
               exit(EXIT_FAILURE);
            }

            /* update index to previous state */
            t_0--;
         } break;

         /* I STATE to {M,I} */
         case I_ST:  /* I connects from M,I at (q_1, t_0) */
         {
            /* current state */
            sc_cur = ISMX(qx0, tx0);
            /* No valid alignment goes to -INF */
            if ( sc_cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible I_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }

            /* possible previous states */
            prv_M = MSMX(qx1, tx0);
            prv_I = ISMX(qx1, tx0);

            /* find maximum next state score */
            if ( sc_max < prv_M ) {
               st_cur = M_ST;
               sc_max = prv_M;
            }
            if ( sc_max < prv_I ) {
               st_cur = I_ST;
               sc_max = prv_I;
            }
            if ( sc_max == -INF ) {
               fprintf( stderr, "ERROR: Failed to trace from I_ST at (%d,%d)\n", t_0, q_0);
               exit(EXIT_FAILURE);
            }

            /* update index to previous state */
            q_0--;
         } break;

         /* N STATE to {N,S} */
         case N_ST:  /* N connects from S,N */
         {
            /* current state */
            sc_cur = XMX(SP_N, q_0);
            /* No valid alignment goes to -INF */
            if ( sc_cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible N_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }

            /* if at beginning of query sequence, then alignment completes at S state, else N state */
            if ( q_0 <= 0 ) {
               st_cur = S_ST;
            } 
            else {
               st_cur = N_ST;
            }
         } break;

         /* B STATE to {N,J} */
         case B_ST:  /* B connects from N, J */
         {
            /* current state */
            sc_cur = XMX(SP_B, q_0);

            /* possible previous states */
            prv_N = XMX(SP_N, q_0) + XSC(SP_N, SP_MOVE);
            prv_J = XMX(SP_J, q_0) + XSC(SP_J, SP_MOVE);

            /* find maximum next state score */            
            if ( sc_max < prv_N ) {
               st_cur = N_ST;
               sc_max = prv_N;
            }
            if ( sc_max < prv_J ) {
               st_cur = J_ST;
               sc_max = prv_J;
            }
            if ( sc_max == -INF ) {
               fprintf( stderr, "ERROR: Failed to trace from B_ST at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }
         } break;

         /* J STATE to {J,E} */
         case J_ST:  /* J connects from E(q_0) or J(q_1) */
         {
            /* current state */
            sc_cur = XMX(SP_J, q_0);
            /* No valid alignment goes to -INF */
            if ( sc_cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible J_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }

            /* possible previous states */
            prv_J = XMX(SP_J, q_1) + XSC(SP_J, SP_LOOP);
            prv_E = XMX(SP_E, q_0) + XSC(SP_E, SP_LOOP);

            /* find maximum next state score */            
            if ( sc_max < prv_J ) {
               st_cur = J_ST;
               sc_max = prv_J;
            }
            if ( sc_max < prv_E ) {
               st_cur = E_ST;
               sc_max = prv_E;
            }
            if ( sc_max == -INF ) {
               fprintf( stderr, "ERROR: Failed to trace from B_ST at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }
         } break;

         default:
         {
            fprintf( stderr, "ERROR: Hit Bogus State!!!\n");
            exit(EXIT_FAILURE);
         }
      }

      fprintf(fp, "[%ld] { %s, (%d, %d) }: %11.4f %11.4f\n", 
         aln->traces->N, STATE_NAMES[st_cur], q_0, t_0, sc_cur, log(sc_cur) );
      ALIGNMENT_Append( aln, st_cur, q_0, t_0 );

      /* For {N,C,J}, we deferred q_0 decrement. */
      if ( (st_cur == N_ST || st_cur == J_ST || st_cur == C_ST) && (st_cur == st_prv) ) {
         q_0--;
      }

      /* Update previous state */
      st_prv = st_cur;
   }

   /* reverse order of traceback */
   ALIGNMENT_Reverse( aln );

   /* find end and begin alignment points (first and last match state) */
   int N  = aln->traces->N;
   for (int i = 0; i < N; i++) {
      if ( tr[i].st == B_ST ) {
         VECTOR_INT_Pushback( aln->tr_beg, i + 1 );
      }
      if ( tr[i].st == E_ST ) {
         VECTOR_INT_Pushback( aln->tr_end, i - 1 );
      }
   }
   aln->beg = VEC_X( aln->tr_beg, 0 );
   aln->end = VEC_X( aln->tr_end, 0 );

   #if DEBUG
   {
      MATRIX_2D* cloud_MX = debugger->cloud_MX;
      MATRIX_2D_Reuse( cloud_MX, Q+1, T+1 );
      MATRIX_2D_Fill( cloud_MX, 0 );
      for ( int i = 0; i < N; i++ ) {
         if ( tr[i].st == M_ST || tr[i].st == I_ST || tr[i].st == D_ST )
            MX_2D( cloud_MX, tr[i].q_0, tr[i].t_0 ) = -1.0;
      }
   }
   #endif

   fclose(fp);

   return STATUS_SUCCESS;
}

/*! FUNCTION:  run_Optimal_Accuracy_Sparse()
 *  SYNOPSIS:  Using <...post> dp matrix, compute optimal accuraccy
 *             matrix which can be used to quickly traceback the 
 *             optimal alignment, stored in <...opt> dp matrix. 
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_Posterior_Optimal_Accuracy_Sparse(    const SEQUENCE*      query,         /* query sequence */
                                          const HMM_PROFILE*   target,        /* target HMM model */
                                          const int            Q,             /* query length */
                                          const int            T,             /* target length */
                                          EDGEBOUNDS*          edg,           /* edgebounds */
                                          RANGE*               dom_range,     /* OPTIONAL: domain range for computing fwd/bck on specific domain. If NULL, computes complete fwd/bck. */
                                          MATRIX_3D_SPARSE*    st_SMX_post,   /* posterior normal state matrix */
                                          MATRIX_2D*           sp_MX_post,    /* posterior special state matrix */
                                          MATRIX_3D_SPARSE*    st_SMX_opt,    /* OUTPUT: optimal normal state matrix */
                                          MATRIX_2D*           sp_MX_opt,     /* OUTPUT: optimal special state matrix */        
                                          float*               sc_final )     /* OUTPUT: final score */
{
   /* vars for matrix access */
   MATRIX_3D_SPARSE*    st_SMX;              /* normal state matrix */
   MATRIX_2D*           sp_MX;               /* special state matrix */

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
   RANGE                T_range;
   RANGE                Q_range;
   bool                 is_q_0_in_dom_range;
   bool                 is_q_1_in_dom_range;

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
      st_SMX      = st_SMX_opt;
      sp_MX       = sp_MX_opt;
      
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

   /* init index */
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
      XMX(SP_B, q_0) = 0.0f; 

      /* only compute if in domain range */
      if ( is_q_0_in_dom_range == true )
      {
         /* FOR every BOUND in zero row (-INF values set during initialization, so unneccessary) */
         for (r_0 = r_0b; r_0 < r_0e; r_0++)
         {
            /* get bound data */
            bnd   = &EDG_X(edg, r_0);
            id    = bnd->id;
            lb_0  = (bnd->lb, 1);           /* can't overflow the left edge */
            rb_0  = (bnd->rb, T);           /* can't overflow the right edge */

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

         /* init lookback 1 row */
         r_1b = r_0b;
         r_1e = r_0e;
      }
   }
   
   /* MAIN RECURSION */
   /* FOR every position in QUERY sequence (row in matrix) */
   for (q_0 = 1; q_0 <= Q; q_0++)
   {
      q_1 = q_0 - 1;

      /* check if query position is in domain */
      is_q_0_in_dom_range = IS_IN_RANGE( Q_range.beg, Q_range.end, q_0 );
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

            /* MAIN RECURSION */
            /* FOR every position in TARGET profile */
            for (t_0 = lb_0; t_0 < rb_0 - 1; t_0++)
            {
               t_1 = t_0 - 1; 
               tx0 = t_0 - bnd->lb;
               tx1 = tx0 - 1;

               /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
               /* best previous state transition (match takes the diag element of each prev state) */
               prv_M =  TSC_DELTA(t_1, M2M) * 
                        (MSMX_X(st_SMX_opt, qx1, tx1) + MSMX_X(st_SMX_post, qx0, tx0));
               prv_I =  TSC_DELTA(t_1, I2M) * 
                        (ISMX_X(st_SMX_opt, qx1, tx1) + MSMX_X(st_SMX_post, qx0, tx0));
               prv_D =  TSC_DELTA(t_1, D2M) * 
                        (DSMX_X(st_SMX_opt, qx1, tx1) + MSMX_X(st_SMX_post, qx0, tx0));
               prv_B =  TSC_DELTA(t_1, B2M) * 
                        ( XMX_X(sp_MX_opt,  B2M, q_1) + MSMX_X(st_SMX_post, qx0, tx0));
               /* best-to-match */
               prv_sum = MAX( MAX( prv_M, prv_I ),
                              MAX( prv_B, prv_D ) );
               MSMX_X(st_SMX_opt, qx0, tx0) = prv_sum;

               /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
               /* previous states (match takes the previous row (upper) of each state) */
               prv_M =  TSC_DELTA(t_0, M2I) * 
                        (MSMX_X(st_SMX_opt, qx1, tx0) + ISMX_X(st_SMX_post, qx0, tx0));
               prv_I =  TSC_DELTA(t_0, I2I) * 
                        (ISMX_X(st_SMX_opt, qx1, tx0) + ISMX_X(st_SMX_post, qx0, tx0));
               /* best-to-insert */
               prv_sum = MAX( prv_M, prv_I );
               ISMX_X(st_SMX_opt, qx0, tx0) = prv_sum;

               /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) */
               /* previous states (match takes the previous column (left) of each state) */
               prv_M =  TSC_DELTA(t_0, M2D) * 
                        (MSMX_X(st_SMX_opt, qx0, tx1) + DSMX_X(st_SMX_post, qx0, tx0));
               prv_D =  TSC_DELTA(t_0, D2D) * 
                        (DSMX_X(st_SMX_opt, qx0, tx1) + DSMX_X(st_SMX_post, qx0, tx0));
               /* best-to-delete */
               prv_sum = MAX( prv_M, prv_D );
               DSMX_X(st_SMX_opt, qx0, tx0) = prv_sum;

               /* UPDATE E STATE */
               prv_M = MSMX(qx0, tx0) + sc_E;
               prv_E = XMX_X(sp_MX_post, SP_E, q_0);
               /* best-to-end */
               prv_sum = MAX( prv_E, prv_M );
               XMX_X(sp_MX_opt, SP_E, q_0) = prv_sum;

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

            /* UNROLLED FINAL LOOP ITERATION */
            {
               t_0 = rb_0 - 1;
               t_1 = t_0 - 1; 
               tx0 = t_0 - bnd->lb;
               tx1 = tx0 - 1;

               /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
               /* best previous state transition (match takes the diag element of each prev state) */
               prv_M =  TSC_DELTA(t_1, M2M) * 
                        (MSMX_X(st_SMX_opt, qx1, tx1) + MSMX_X(st_SMX_post, qx0, tx0));
               prv_I =  TSC_DELTA(t_1, I2M) * 
                        (ISMX_X(st_SMX_opt, qx1, tx1) + MSMX_X(st_SMX_post, qx0, tx0));
               prv_D =  TSC_DELTA(t_1, D2M) * 
                        (DSMX_X(st_SMX_opt, qx1, tx1) + MSMX_X(st_SMX_post, qx0, tx0));
               prv_B =  TSC_DELTA(t_1, B2M) * 
                        ( XMX_X(sp_MX_opt,  B2M, q_1) + MSMX_X(st_SMX_post, qx0, tx0));
               prv_B = XMX(SP_B, q_1) + TSC(t_1, B2M); /* from begin match state (new alignment) */
               /* best-to-match */
               prv_sum = MAX( MAX( prv_M, prv_I ),
                              MAX( prv_B, prv_D ) );
               MSMX_X(st_SMX_opt, qx0, tx0) = prv_sum;

               /* FIND SUM OF PATHS TO INSERT STATE (unrolled) */
               ISMX(qx0, tx0) = -INF;

               /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) */
               /* previous states (match takes the previous column (left) of each state) */
               prv_M =  TSC_DELTA(t_0, M2D) * 
                        (MSMX_X(st_SMX_opt, qx0, tx1) + DSMX_X(st_SMX_post, qx0, tx0));
               prv_M =  TSC_DELTA(t_0, D2D) * 
                        (DSMX_X(st_SMX_opt, qx0, tx1) + DSMX_X(st_SMX_post, qx0, tx0));
               /* best-to-delete */
               prv_sum = MAX( prv_M, prv_D );
               DSMX_X(st_SMX_opt, qx0, tx0) = prv_sum;

               /* UPDATE E STATE (unrolled) */
               prv_M = MSMX_X(st_SMX_opt, qx0, tx0);
               prv_D = DSMX_X(st_SMX_opt, qx0, tx0);
               prv_E = XMX_X(sp_MX_opt, SP_E, q_0);
               /* best-to-end */
               prv_sum = MAX( MAX(  prv_D, prv_M ),
                                    prv_E ); 
               XMX_X(sp_MX_opt, SP_E, q_0) = prv_sum;

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
      /* J state ( J->J or E->J (E's loop) ) */
      prv_J =  XSC_DELTA(SP_J, SP_LOOP) * 
               (XMX_X(sp_MX_opt, SP_J, q_1) + XMX_X(sp_MX_post, SP_J, q_0));
      prv_E =  XSC_DELTA(SP_E,SP_LOOP) * 
               (XMX_X(sp_MX_opt, SP_E, q_0));
      prv_sum = MAX( prv_J, prv_E );
      XMX_X(sp_MX_opt, SP_J, q_0) = prv_sum;   

      /* C state ( C->C or E->C ) */
      prv_C =  XSC_DELTA(SP_C, SP_LOOP) *
               (XMX_X(sp_MX_opt, SP_C, q_1) + XMX_X(sp_MX_post, SP_C, q_0));
      prv_E =  XSC_DELTA(SP_E,SP_MOVE) *
               (XMX_X(sp_MX_opt, SP_E, q_0));
      prv_sum = MAX( prv_C, prv_E );
      XMX(SP_C, q_0) = prv_sum;

      /* N state ( N->N (tail) ) */
      prv_N =  XSC_DELTA(SP_N, SP_LOOP) *
               (XMX_X(sp_MX_opt, SP_N, q_1) + XMX_X(sp_MX_post, SP_N, q_0));
      XMX(SP_N, q_0) = prv_N;

      /* B state ( N->B (N's move) or J->B (J's move) ) */
      prv_N =  XSC_DELTA(SP_N, SP_MOVE) *
               (XMX_X(sp_MX_opt, SP_N, q_0));
      prv_J =  XSC_DELTA(SP_N, SP_MOVE) *
               (XMX_X(sp_MX_opt, SP_J, q_0));
      prv_sum = MAX( prv_N, prv_J ); 
      XMX(SP_B, q_0) = prv_sum;   

      /* SET CURRENT ROW TO PREVIOUS ROW */
      r_1b = r_0b;
      r_1e = r_0e;
   }

   /* T state */
   sc_best     = XMX_X(sp_MX_post, SP_C, Q);
   *sc_final   = sc_best;

   return STATUS_SUCCESS;
}