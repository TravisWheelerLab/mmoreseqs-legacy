/*******************************************************************************
 *  - FILE:    	viterbi_linear.h
 *  - DESC:  	The Viterbi Algorithm and Traceback for Sequence Alignment Search.
 *******************************************************************************/

#ifndef _VITERBI_LINEAR_H
#define _VITERBI_LINEAR_H

/*! FUNCTION:  run_Viterbi_Linear()
 *  SYNOPSIS:  Run Viterbi Algorithm (Seq-to-Profile, general unoptimized)
 *  ARGS:      <query>     query sequence,
 *             <target>    HMM model,
 *             <Q>         Query length,
 *             <T>         Target length,
 *             <st_MX>     Normal state matrix,
 *             <sp_MX>     Special state matrix
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Viterbi_Linear(const SEQUENCE* query,
                       const HMM_PROFILE* target,
                       const int Q,
                       const int T,
                       MATRIX_3D* st_MX,
                       MATRIX_2D* sp_MX,
                       float* sc_final);

#endif /* _VITERBI_LINEAR_H */
