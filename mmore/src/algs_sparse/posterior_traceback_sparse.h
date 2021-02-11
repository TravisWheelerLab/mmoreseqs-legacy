/*******************************************************************************
 *     FILE:   posterior_traceback_sparse.h
 *  PURPOSE:   The Posterior Probability and Optimal Alignment.
 *
 *   AUTHOR:   Dave Rich
 *******************************************************************************/

#ifndef _POSTERIOR_TRACEBACK_SPARSE_H
#define _POSTERIOR_TRACEBACK_SPARSE_H

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
                                            MATRIX_3D_SPARSE*       st_SMX_post,      /* posterior normal matrix */
                                            MATRIX_2D*              sp_MX_post,       /* posterior special matrix */
                                            MATRIX_3D_SPARSE*       st_SMX_opt,       /* optimal accuracy normal matrix */
                                            MATRIX_2D*              sp_MX_opt,        /* optimal accuracy special matrix */
                                            ALIGNMENT*              aln );            /* OUTPUT: optimal alignment */

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
                                          RANGE*               in_dom_range,  /* OPTIONAL: domain range for computing fwd/bck on specific domain. If NULL, computes complete fwd/bck. */
                                          MATRIX_3D_SPARSE*    st_SMX_post,   /* posterior normal state matrix */
                                          MATRIX_2D*           sp_MX_post,    /* posterior special state matrix */
                                          MATRIX_3D_SPARSE*    st_SMX_opt,    /* OUTPUT: optimal normal state matrix */
                                          MATRIX_2D*           sp_MX_opt,     /* OUTPUT: optimal special state matrix */        
                                          float*               sc_final );    /* OUTPUT: final score */

#endif /* _POSTERIOR_TRACEBACK_SPARSE_H */
