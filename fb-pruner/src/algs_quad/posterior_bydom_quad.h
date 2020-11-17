/*******************************************************************************
 *     FILE:   posterior_bydom_quad.h
 *  PURPOSE:   The Posterior Probability and Optimal Alignment.
 *
 *   AUTHOR:   Dave Rich
 *******************************************************************************/

#ifndef _POSTERIOR_BYDOM_QUAD_H
#define _POSTERIOR_BYDOM_QUAD_H


/*! FUNCTION:  run_Posterior_ByDomain_Quad()
 *  SYNOPSIS:  Filled dp matrices for forward <st_MX_fwd> and backward <st_MX_bck>.
 *             Compute the Posterior Probability by multiplying probabilities (added in log space) of Forward and Backward.
 *             Results stored in supplied <st_MX_post> (can override input matrices).
 *             NOTE: Modeled after <p7_Decoding()> and <p7_DomainDecoding()>.
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int 
run_Posterior_ByDomain_Quad(  SEQUENCE*      q_seq,            /* query sequence */
                              HMM_PROFILE*   t_prof,           /* target hmm model */
                              int            Q,                /* query length */
                              int            T,                /* target length */
                              HMM_BG*        bg,               /* hmm background model */
                              MATRIX_3D*     st_MX_fwd,        /* normal state matrix for forward */
                              MATRIX_2D*     sp_MX_fwd,        /* special state matrix for forward */
                              MATRIX_3D*     st_MX_bck,        /* normal state matrix for backward */
                              MATRIX_2D*     sp_MX_bck,        /* special state matrix for backward */
                              MATRIX_3D*     st_MX_post,       /* OUTPUT: normal state matrix for posterior */
                              MATRIX_2D*     sp_MX_post,       /* OUTPUT: special state matrix for posterior */     
                              DOMAIN_DEF*    dom_def );        /* OUTPUT: domain data */


/*! FUNCTION:  run_Decode_Normal_Posterior_Quad()
 *  SYNOPSIS:  Using <...fwd> and <...bck> dp matrices to create special state posterior into <...post>.
 *             Can store matrix in <...fwd> or <...bck>.
 *             NOTE: Modeled after <p7_Decoding()> and <p7_DomainDecoding()>
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int
run_Decode_Posterior_ByDomain_Quad( SEQUENCE*         q_seq,            /* query sequence */
                                    HMM_PROFILE*      t_prof,           /* target hmm model */
                                    int               Q,                /* full query size */
                                    int               T,                /* full target size */
                                    RANGE*            Q_range,          /* query range */
                                    RANGE*            T_range,          /* target range */
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
run_Decode_Special_Posterior_ByDomain_Quad(     SEQUENCE*         q_seq,            /* query sequence */
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
run_Rescore_Isolated_ByDomain(  SEQUENCE*         q_seq,            /* query sequence */
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
run_Null2_ByExpectation_ByDomain(   SEQUENCE*         query,            /* query sequence */
                                    HMM_PROFILE*      target,           /* target hmm model */
                                    int               Q,                /* query length */
                                    int               T,                /* target length */
                                    RANGE*            Q_range,          /* query range */
                                    RANGE*            T_range,          /* target range */ 
                                    MATRIX_3D*        st_MX_post,       /* posterior normal matrix */
                                    MATRIX_2D*        sp_MX_post,       /* posterior special matrix */
                                    DOMAIN_DEF*       dom_def );        /* OUTPUT: domain def's null2_sc vector */


#endif /* _POSTERIOR_QUAD_H */
