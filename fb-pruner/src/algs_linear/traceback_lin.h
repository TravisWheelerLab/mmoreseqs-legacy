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

/* === FUNCTIONS === */
int run_Traceback_Quad(    const SEQUENCE*     query,
                           const HMM_PROFILE*  target,
                           const int           Q, 
                           const int           T,
                           MATRIX_3D*          st_MX,
                           MATRIX_2D*          sp_MX,
                           ALIGNMENT*          aln );

int traceback_Append(   ALIGNMENT*  	aln,
					         TRACE* 		   tr,
                        const int   	st,
                        const int   	i,
                        const int   	j );

int traceback_Reverse(  ALIGNMENT* aln );

#endif /* _VITERBI_H */