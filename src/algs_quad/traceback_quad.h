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

/* === INCLUDES === */
// #include "objects/sequence.h"
// #include "objects/hmm_profile.h"
// #include "objects/alignment.h"

/* === FUNCTIONS === */
void traceback_Build(const SEQUENCE*     query,
                     const HMM_PROFILE*  target,
                     const int           Q, 
                     const int           T,
                     float*              st_MX,
                     float*              sp_MX,
                     ALIGNMENT*          aln);

void traceback_Append(ALIGNMENT*  aln,
                      int          st,
                      int          i,
                      int          j);

void traceback_Reverse(ALIGNMENT* aln);

#endif /* _VITERBI_H */