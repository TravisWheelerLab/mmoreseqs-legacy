/*******************************************************************************
 *     FILE:   posterior_quad.h
 *  PURPOSE:   The Posterior Probability and Optimal Alignment.
 *
 *   AUTHOR:   Dave Rich
 *******************************************************************************/

#ifndef _POSTERIOR_QUAD_H
#define _POSTERIOR_QUAD_H

/*! FUNCTION:  run_MaxPost_Quad()
 *  SYNOPSIS:  Compute the Max Posterior Probability and obtain optimal alignment.
 *             Computes the Forward and Backward.  Max Posterior computed by Viterbi.
 *
 *  RETURN:    <STATUS_SUCCESS> if no errors.
 */
int 
run_MaxPost_Quad( SEQUENCE*         query,            /* query sequence */
                  HMM_PROFILE*      target,           /* target HMM model */
                  int               Q,                /* query length */
                  int               T,                /* target length */
                  MATRIX_3D*        st_MX_fwd,        /* normal matrix for forward */
                  MATRIX_2D*        sp_MX_fwd,        /* special matrix for forward */
                  MATRIX_3D*        st_MX_bck,        /* normal matrix for backward */
                  MATRIX_2D*        sp_MX_bck,        /* special matrix for backward */
                  MATRIX_3D*        st_MX_post,       /* OUTPUT: normal matrix for posterior (can overwrite fwd and bck data) */
                  MATRIX_2D*        sp_MX_post,       /* OUTPUT: special matrix for posterior (can overwrite fwd and bck data) */
                  ALIGNMENT*        aln,              /* OUTPUT: alignment */
                  DOMAIN_DEF*       dom_def );        /* OUTPUT: domain definition scores */


/*! FUNCTION:  run_Posterior_Quad()
 *  SYNOPSIS:  Filled dp matrices for forward <st_MX_fwd> and backward <st_MX_bck>.
 *             Compute the Posterior Probability by multiplying probabilities (added in log space) of Forward and Backward.
 *             Results stored in supplied <st_MX_post> (can override input matrices).
 *             NOTE: Modeled after <p7_Decoding()> and <p7_DomainDecoding()>.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int 
run_Posterior_Quad(  SEQUENCE*      q_seq,            /* query sequence */
                     HMM_PROFILE*   t_prof,           /* target hmm model */
                     int            Q,                /* query length */
                     int            T,                /* target length */
                     MATRIX_3D*     st_MX_fwd,        /* normal state matrix for forward */
                     MATRIX_2D*     sp_MX_fwd,        /* special state matrix for forward */
                     MATRIX_3D*     st_MX_bck,        /* normal state matrix for backward */
                     MATRIX_2D*     sp_MX_bck,        /* special state matrix for backward */
                     MATRIX_3D*     st_MX_post,       /* OUTPUT: normal state matrix for posterior */
                     MATRIX_2D*     sp_MX_post,       /* OUTPUT: special state matrix for posterior */     
                     DOMAIN_DEF*    dom_def );        /* OUTPUT: domain data */


/*! FUNCTION:  test_Multidomain_Region()
 *  SYNOPSIS:  Returns whether sequence range <q_beg,q_end> contains multiple domains.
 *             NOTES: More precisely: return TRUE if  \max_z [ \min (B(z), E(z)) ]  >= rt3
 *             where
 *             E(z) = expected number of E states occurring in region before z is emitted
 *                = \sum_{y=i}^{z} eocc[i]  =  etot[z] - etot[i-1]
 *             B(z) = expected number of B states occurring in region after z is emitted
 *                = \sum_{y=z}^{j} bocc[i]  =  btot[j] - btot[z-1]
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
bool
test_Multidomain_Region(   DOMAIN_DEF*    dom_def,
                           int            q_beg,
                           int            q_end );


/*! FUNCTION:  run_Decode_Normal_Posterior_Quad()
 *  SYNOPSIS:  Using <...fwd> and <...bck> dp matrices to create special state posterior into <...post>.
 *             Can store matrix in <...fwd> or <...bck>.
 *             NOTE: Modeled after <p7_Decoding()> and <p7_DomainDecoding()>
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_Decode_Posterior_Quad( SEQUENCE*         q_seq,            /* query sequence */
                           HMM_PROFILE*      t_prof,           /* target hmm model */
                           int               Q,                /* query length */
                           int               T,                /* target length */
                           MATRIX_3D*        st_MX_fwd,        /* normal state matrix for forward */
                           MATRIX_2D*        sp_MX_fwd,        /* special state matrix for forward */
                           MATRIX_3D*        st_MX_bck,        /* normal state matrix for backward */
                           MATRIX_2D*        sp_MX_bck,        /* special state matrix for backward */
                           MATRIX_3D*        st_MX_post,       /* OUTPUT: normal state matrix for posterior */
                           MATRIX_2D*        sp_MX_post );     /* OUTPUT: normal state matrix for posterior */


/*! FUNCTION:  run_Decode_Posterior_Quad()
 *  SYNOPSIS:  Using <...fwd> and <...bck> dp matrices to create special state posterior into <...post>.
 *             Can store matrix in <...fwd> or <...bck>.
 *             NOTE: Modeled after <p7_Decoding()> and <p7_DomainDecoding()>
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_Decode_Special_Posterior_Quad(  SEQUENCE*         q_seq,            /* query sequence */
                                    HMM_PROFILE*      t_prof,           /* target hmm model */
                                    int               Q,                /* query length */
                                    int               T,                /* target length */
                                    MATRIX_2D*        sp_MX_fwd,        /* special state matrix for forward */
                                    MATRIX_2D*        sp_MX_bck,        /* special state matrix for backward */
                                    MATRIX_2D*        sp_MX_post );     /* OUTPUT: special state matrix for posterior */


/*! FUNCTION:  run_Decode_Posterior_Quad()
 *  SYNOPSIS:  Using posterior special state matrix <sp_MX_post> to compute domains.
 *             NOTE: Modeled after <p7_DomainDecoding()>
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_DecodeDomain_Posterior_Quad( SEQUENCE*         q_seq,            /* query sequence */
                                 HMM_PROFILE*      t_prof,           /* target hmm model */
                                 int               Q,                /* query length */
                                 int               T,                /* target length */
                                 MATRIX_2D*        sp_MX_post,       /* special state matrix for posterior */ 
                                 DOMAIN_DEF*       dom_def );        /* OUTPUT: domain data */

/*! FUNCTION:  run_DecodeDomain_Quad()
 *  SYNOPSIS:  Using posterior special state matrix <sp_MX_bck> and <sp_MX_fwd> to compute domains.
 *             NOTE: Modeled after <p7_DomainDecoding()>
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_DecodeDomain_Quad(  SEQUENCE*         q_seq,            /* query sequence */
                        HMM_PROFILE*      t_prof,           /* target hmm model */
                        int               Q,                /* query length */
                        int               T,                /* target length */
                        MATRIX_2D*        sp_MX_fwd,       /* special state matrix for forward */ 
                        MATRIX_2D*        sp_MX_bck,       /* special state matrix for backward */ 
                        DOMAIN_DEF*       dom_def );        /* OUTPUT: domain data */  

/*! FUNCTION:  run_Rescore_Isolated_Domain()
 *  SYNOPSIS:  
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int 
run_Rescore_Isolated_Domain(  SEQUENCE*         q_seq,            /* query sequence */
                              HMM_PROFILE*      t_prof,           /* target hmm model */
                              RANGE*            Q_range,          /* query range */
                              int               T,                /* target length */
                              MATRIX_3D*        st_MX_fwd,        /* normal state matrix for forward */
                              MATRIX_2D*        sp_MX_fwd,        /* special state matrix for forward */
                              MATRIX_3D*        st_MX_bck,        /* normal state matrix for backward */
                              MATRIX_2D*        sp_MX_bck,        /* special state matrix for backward */
                              ALIGNMENT*        aln,              /* OUTPUT: domain alignment */
                              DOMAIN_DEF*       dom_def );        /* OUTPUT: domain data */


/*! FUNCTION:  run_Null2_By_Expectation()
 *  SYNOPSIS:  Modeled after HMMER p7_GNull2_ByExpectation().
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_Null2_ByExpectation(   SEQUENCE*         query,            /* query sequence */
                           HMM_PROFILE*      target,           /* target hmm model */
                           RANGE*            Q,                /* query length */
                           int               T,                /* target length */
                           MATRIX_3D*        st_MX_post,       /* posterior normal matrix */
                           MATRIX_2D*        sp_MX_post,       /* posterior special matrix */
                           DOMAIN_DEF*       dom_def );        /* OUTPUT: domain def's null2_sc vector */


// /*! FUNCTION:   run_Posterior_Traceback_Quad()
//  *  SYNOPSIS:   Run Viterbi Traceback to recover Optimal Alignment.
//  *
//  *    RETURN:   Return <STATUS_SUCCESS> if no errors.
//  */
// int 
// run_Posterior_Traceback_Quad(    const SEQUENCE*     query,       /* query sequence */
//                                  const HMM_PROFILE*  target,      /* HMM model */
//                                  const int           Q,           /* query/seq length */
//                                  const int           T,           /* target/model length */
//                                  MATRIX_3D*          st_MX,       /* Normal State (Match, Insert, Delete) Matrix */
//                                  MATRIX_2D*          sp_MX,       /* Special State (J,N,B,C,E) Matrix */
//                                  ALIGNMENT*          aln );       /* OUTPUT: Traceback Alignment */


#endif /* _POSTERIOR_QUAD_H */
