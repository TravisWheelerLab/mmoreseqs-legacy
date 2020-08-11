/*******************************************************************************
 *
 *  FILE:    traceback.h
 *  PURPOSE: Traceback for Viterbi Algorithm.
 *
 *  AUTHOR:  Dave Rich
 *
 *******************************************************************************/

#ifndef _VITERBI_H
#define _VITERBI_H

/*
 *  FUNCTION:  run_Traceback_Quad()
 *  SYNOPSIS:  Run Viterbi Traceback to recover Optimal Alignment.
 *
 *  ARGS:      <query>     query sequence,
 *             <target>    HMM model,
 *             <Q>         query/seq length,
 *             <T>         target/model length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Maalnix,
 *             <sp_MX>     Special State (J,N,B,C,E) Maalnix,
 *             <aln>       Traceback Alignment
 *
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int run_Traceback_Lin(  const SEQUENCE*     query,
                        const HMM_PROFILE*  target,
                        const int           Q, 
                        const int           T,
                        MATRIX_3D*          st_MX,
                        MATRIX_2D*          sp_MX,
                        ALIGNMENT*          aln  );


#endif /* _VITERBI_H */