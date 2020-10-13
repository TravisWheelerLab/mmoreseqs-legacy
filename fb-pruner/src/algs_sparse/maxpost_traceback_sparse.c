/*******************************************************************************
 *    FILE:       traceback_sparse.c
 *    PURPOSE:    Traceback for Viterbi Algorithm for Sparse Matrices.
 *
 *    AUTHOR:     Dave Rich
 *    BUG:     
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
#include "../algs_quad/algs_quad.h"


/* header */
#include "maxpost_traceback_sparse.h"

/*  FUNCTION:  run_MaxExp_Traceback_Sparse_2()
 *  SYNOPSIS:  Run Viterbi Traceback to recover Optimal Alignment.
 *             Version 2: My implementation. Verifies that Alignment agrees with Matrix data.
 *
 *    RETURN:  Return <STATUS_SUCCESS> if no errors.
 */
int run_MaxPost_Traceback_Sparse(   const SEQUENCE*      query,      /* query sequence */
                                    const HMM_PROFILE*   target,     /* HMM model */
                                    const int            Q,          /* query/seq length */
                                    const int            T,          /* target/model length */
                                    MATRIX_3D_SPARSE*    st_SMX,     /* Normal State (Match, Insert, Delete) Matrix */
                                    MATRIX_2D*           sp_MX,      /* Special State (J,N,B,C,E) Matrix */
                                    EDGEBOUNDS*          edg,        /* edgebounds of sparse matrix */
                                    ALIGNMENT*           aln )       /* OUTPUT: Traceback Alignment */
{
   /* vars for accessing query/target data structs */
   char     a;                               /* store current character in sequence */
   int      A;                               /* store int value of character */
   char*    seq;                             /* alias for getting seq */
   bool     is_local;                        /* whether */

   /* vars for indexing into data matrices by row-col */
   int      b, d, i, j, k;                   /* antidiagonal, row, column indices */
   int      q_0, q_1;                        /* real index of current and previous rows (query) */
   int      qx0, qx1;                        /* maps column index into data index (query) */
   int      t_0, t_1;                        /* real index of current and previous columns (target) */
   int      tx0, tx1;                        /* maps target index into data index (target)  */
   int      t_range;                         /* range of targets on current row */

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
   bool     rb_T;                            /* checks if edge touches right bound of matrix */

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
   int            st_cur, st_prv;            /* current, previous state in traceback */
   int            q_prv, t_prv;              /* previous  */
   float          max_sc;                    /* maximum score */
   int            max_st;                    /* maximum state */
   TRACE*         tr;                        /* trace object for appending */
   const float    tol = 1e-5;                /* acceptable tolerance range for "equality tests" */

   /* --------------------------------------------------------------------------- */

   /* get sequence and trace */
   seq   = query->seq; 
   tr    = aln->traces->data; 

   /* intial cell */
   q_0 = q_prv = Q;
   t_0 = t_prv = 0;

   q_prv = -1;
   t_prv = -1;

   /* add every edgebound from current row */
   r_0b = edg->N - 1;
   r_0 = r_0e = r_0b;
   while ( r_0 >= 0 && EDG_X(edg, r_0).id >= q_0 ) {
      r_0--;
   }
   r_0e = r_0;

   /* local or global? */
   is_local = target->isLocal;

   /* clear memory for trace */
   ALIGNMENT_Reuse( aln, Q, T );

   /* Backtracing, so C is end state */
   ALIGNMENT_Append( aln, tr, T_ST, q_0, t_0 );
   ALIGNMENT_Append( aln, tr, C_ST, q_0, t_0 );
   st_prv = C_ST;

   /* End of trace is S state */
   while (st_prv != S_ST)
   {
      /* if we have reached the end of the query sequence */
      if (q_0 == 0) {
         ALIGNMENT_Append( aln, tr, S_ST, q_0, t_0 );
         break;
      } 

      /* if q_0 has been decremented, then get edgebound range of next row  */
      if ( q_prv != q_0 ) {
         r_0b = r_0 = r_0e;
         while ( r_0 >= 0 && EDG_X(edg, r_0).id >= q_0 ) {
            r_0--;
         }
         r_0e = r_0;
         /* reset r_0 to beginning of row */
         r_0 = r_0b;
      }

      printf("r_0=(%d,%d)\n", r_0b, r_0e);

      /* TODO: can be optimized by only updating during certain cases */
      /* if either t_0 or q_0 have been decremented, then lookup index offset */
      if ( t_prv != t_0 || q_prv != q_0 ) {
         /* search through bounds in row for bound containing index */
         for ( r_0 = r_0b; r_0 > r_0e; r_0-- ) {
            /* find bound that contains current index */
            bnd = &EDG_X(edg, r_0);
            if ( t_0 >= bnd->lb && t_0 < bnd->rb ) {
               /* fetch data mapping bound start location to data block in sparse matrix */
               qx0 = VECTOR_INT_Get( st_SMX->imap_cur, r_0 );    /* (q_0, t_0) location offset */
               qx1 = VECTOR_INT_Get( st_SMX->imap_prv, r_0 );    /* (q_1, t_0) location offset */

               tx0 = t_0 - bnd->lb;    /* total_offset = offset_location - starting_location */
               tx1 = tx0 - 1;

               break;
            }
         }
      }

      printf("q_0,t_0=(%d,%d)\n", q_0, t_0);
      printf("qx0=%d, qx1=%d\n", qx0, qx1);

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
            cur   = XMX(SP_C, q_0);
            /* max current state */
            max_sc   = -INF;

            if ( cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible C_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }
            
            /* possible previous states */
            prv_C = XMX(SP_C, q_1);
            prv_E = XMX(SP_E, q_0);

            if ( max_sc < prv_C ) {
               st_cur = C_ST;
               max_sc = prv_C;
            }
            if ( max_sc < prv_E ) {
               st_cur = E_ST;
               max_sc = prv_E;
            }

         } break;

         /* E STATE to {M, D} */
         case E_ST:  /* E connects from any M state. t_0 is set here. */
         {
            /* current state */
            cur = XMX(SP_E, q_0);
            /* set if correct match state found */
            int found = false;

            if ( cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible E_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }

            if ( is_local )  /* local mode: ends in M */
            {
               /* can't come from D, in a *local* Viterbi alignment. */
               st_cur = M_ST;
               max_sc = -INF;
               max_st = -1;

               /* possible previous states (any M state) */
               /* FOR every BOUND in current ROW */
               for ( r_0 = r_0b; r_0 > r_0e; r_0-- )
               {
                  /* get bound data */
                  bnd   = &EDG_X(edg, r_0);
                  id    = bnd->id;
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

                     /* verifies if scores agree with true previous state in alignment */
                     if ( max_sc < prv_M ) {
                        max_st = t_0;
                     }
                  }
               }
            }
            else     /* glocal mode: we either come from D_M or M_M */
            {
               /* can't come from D, in a *local* Viterbi alignment. */
               st_cur = M_ST;
               max_sc = -INF;
               max_st = -1;

               /* get bound data */
               bnd   = &EDG_X(edg, r_0b);
               id    = bnd->id;
               rb_T  = (bnd->rb > T);
               rb_0  = MIN(bnd->rb, T);         /* can't overflow the right edge */

               /* fetch data mapping bound start location to data block in sparse matrix */
               qx0   = VECTOR_INT_Get( st_SMX->imap_cur, r_0b );    /* (q_0, t_0) location offset */
               tx0   = T - bnd->lb;

               /* possible previous states */
               if ( rb_T ) {
                  prv_M = MSMX(qx0, tx0);
                  prv_D = DSMX(qx0, tx0);
               }
               else {
                  prv_M = -INF;
                  prv_D = -INF;
               }

               /* verifies if scores agree with true previous state in alignment */
               if ( max_sc < prv_M ) {
                  max_sc = prv_M;
                  st_cur = M_ST;
                  t_0 = T;
               }
               if ( max_sc < prv_D ) {
                  max_sc = prv_D;
                  st_cur = D_ST;
                  t_0 = T;
               }
            }
         } break;

         /* M STATE to {B,M,I,D} */
         case M_ST:  /* M connects from (q_1, t_1), or B */
         {
            max_sc = -INF;
            max_st = -1;

            /* current state */
            cur = MSMX(qx0, tx0);

            /* No valid alignment goes to -INF */
            if ( cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible M_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }

            /* possible previous states */
            prv_B = XMX(SP_B, q_1) + TSC(t_1, B2M) + MSC(t_0, A);
            prv_M = MSMX(qx1, tx1) + TSC(t_1, M2M) + MSC(t_0, A);
            prv_I = ISMX(qx1, tx1) + TSC(t_1, I2M) + MSC(t_0, A);
            prv_D = DSMX(qx1, tx1) + TSC(t_1, D2M) + MSC(t_0, A);

            /* verifies if scores agree with true previous state in alignment */
            if ( max_sc < prv_B ) {
               max_sc = prv_B;
               st_cur = B_ST;
            }
            if ( max_sc < prv_M ) {
               max_sc = prv_M;
               st_cur = M_ST;
            }
            if ( max_sc < prv_I ) {
               max_sc = prv_I;
               st_cur = I_ST;
            }
            if ( max_sc < prv_D ) {
               max_sc = prv_D;
               st_cur = D_ST;
            }

            /* update index to previous state */
            t_0--; q_0--;
         } break;

         /* D STATE to {M,D} */
         case D_ST:  /* D connects from M,D at (q_0, t_1) */
         {
            max_sc = -INF;
            max_st = -1;

            /* current state */
            cur = DSMX(qx0, tx0);

            t_1 = t_0 - 1;

            /* No valid alignment goes to -INF */
            if ( cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible D_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }

            /* possible previous states */
            prv_M = MSMX(qx0, tx1) + TSC(t_1, M2D);
            prv_D = DSMX(qx0, tx1) + TSC(t_1, D2D);

            /* verifies if scores agree with true previous state in alignment */
            if ( max_sc < prv_M ) {
               max_sc = prv_M;
               st_cur = M_ST;
            }
            if ( max_sc < prv_D ) {
               max_sc = prv_D;
               st_cur = D_ST;
            }

            /* update index to previous state */
            t_0--;
         } break;

         /* I STATE to {M,I} */
         case I_ST:  /* I connects from M,I at (q_1, t_0) */
         {
            max_sc = -INF;
            max_st = -1;

            /* current state */
            cur = ISMX(qx0, tx0);

            /* No valid alignment goes to -INF */
            if ( cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible I_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }

            /* possible previous states */
            prv_M = MSMX(qx1, tx0) + TSC(t_0, M2I) + ISC(t_0, A);
            prv_I = ISMX(qx1, tx0) + TSC(t_0, I2I) + ISC(t_0, A);

            /* verifies if scores agree with true previous state in alignment */
            if ( max_sc < prv_M ) {
               max_sc = prv_M;
               st_cur = M_ST;
            }
            if ( max_sc < prv_I ) {
               max_sc = prv_I;
               st_cur = I_ST;
            }

            /* update index to previous state */
            q_0--;
         } break;

         /* N STATE to {N,S} */
         case N_ST:  /* N connects from S,N */
         {
            /* current state */
            cur = XMX(SP_N, q_0);

            /* No valid alignment goes to -INF */
            if ( cur == -INF ) {
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
            max_sc = -INF;
            max_st = -1;

            /* current state */
            cur = XMX(SP_B, q_0);

            /* possible previous states */
            prv_N = XMX(SP_N, q_0) + XSC(SP_N, SP_MOVE);
            prv_J = XMX(SP_J, q_0) + XSC(SP_J, SP_MOVE);

            /* verifies that scores agree with true previous state in alignment */
            if ( max_sc < prv_N ) {
               max_sc = prv_N;
               st_cur = N_ST;
            }
            if ( max_sc < prv_J ) {
               max_sc = prv_J;
               st_cur = J_ST;
            }
         } break;

         /* J STATE to {J,E} */
         case J_ST:  /* J connects from E(q_0) or J(q_1) */
         {
            max_sc = -INF;
            max_st = -1;

            /* current state */
            cur = XMX(SP_J, q_0);

            if ( cur == -INF ) {
               fprintf( stderr, "ERROR: Impossible J_ST reached at (%d,%d)\n", q_0, t_0);
               exit(EXIT_FAILURE);
            }

            /* possible previous states */
            prv_J = XMX(SP_J, q_1) + XSC(SP_J, SP_LOOP);
            prv_E = XMX(SP_E, q_0) + XSC(SP_E, SP_LOOP);

            /* verifies that scores agree with true previous state in alignment */
            if ( max_sc < prv_J ) {
               max_sc = prv_J;
               st_cur = J_ST;
            }
            if ( max_sc < prv_E ) {
               max_sc = prv_E;
               st_cur = E_ST;
            }
         } break;

         default:
         {
            fprintf( stderr, "ERROR: Hit Bogus State!!!\n");
            exit(EXIT_FAILURE);
         }
      }

      ALIGNMENT_Append( aln, tr, st_cur, q_0, t_0 );

      /* For {N,C,J}, we deferred q_0 decrement. */
      if ( (st_cur == N_ST || st_cur == J_ST || st_cur == C_ST) && (st_cur == st_prv) ) {
         q_0--;
      }

      /* Update previous state */
      st_prv = st_cur;
   }

   /* reverse order of alnaceback */
   ALIGNMENT_Reverse( aln );

   /* find end and begin alignment points (first and last match state) */
   int N  = aln->traces->N;
   for (int i = 0; i < N; ++i) {
      if ( tr[i].st == B_ST ) {
         VECTOR_INT_Pushback( aln->tr_beg, i + 1 );
      }
      if ( tr[i].st == E_ST ) {
         VECTOR_INT_Pushback( aln->tr_end, i - 1 );
      }
   }
   aln->beg = aln->tr_beg->data[0];
   aln->end = aln->tr_end->data[0];

   #if DEBUG
   {
      MATRIX_2D* cloud_MX = debugger->cloud_MX;
      MATRIX_2D_Reuse( cloud_MX, Q+1, T+1 );
      MATRIX_2D_Fill( cloud_MX, 0 );
      for ( int i = 0; i < N; i++ ) {
         if ( tr[i].st == M_ST || tr[i].st == I_ST || tr[i].st == D_ST )
            MX_2D( cloud_MX, tr[i].i, tr[i].j ) = -1.0;
      }
   }
   #endif

   return STATUS_SUCCESS;
}

