/*******************************************************************************
 *  FILE:      drich_funcs.h
 *  PURPOSE:   David Rich Functions
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef DRICH_EDIT_H
#define DRICH_EDIT_H

#include "easel.h"
#include "esl_alphabet.h" /* ESL_DSQ, ESL_ALPHABET */
#include "esl_dmatrix.h"  /* ESL_DMATRIX           */
#include "esl_getopts.h"  /* ESL_GETOPTS           */
#include "esl_histogram.h"      /* ESL_HISTOGRAM         */
#include "esl_hmm.h"          /* ESL_HMM               */
#include "esl_keyhash.h"        /* ESL_KEYHASH           */
#include "esl_mixdchlet.h"  /* ESL_MIXDCHLET         */
#include "esl_msa.h"    /* ESL_MSA               */
#include "esl_random.h"   /* ESL_RANDOMNESS        */
#include "esl_rand64.h" /* ESL_RAND64 */
#include "esl_sq.h"   /* ESL_SQ                */
#include "esl_scorematrix.h"    /* ESL_SCOREMATRIX       */
#include "esl_stopwatch.h"      /* ESL_STOPWATCH         */

/** FUNCTION:  dp_matrix_Save()
 *  SYNOPSIS:  Save dynamic programming matrix to file.
 *
 *  ARGS:      <Q>         query length,
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix
 *             <f>         Filename
 */
void DP_MATRIX_Dump(  const int   Q, 
                      const int   T,
                      ESL_DSQ*    dsq,
                      P7_PROFILE* gm,
                      P7_GMX*     gx,
                      FILE*       fp );

/** FUNCTION: 
 *
 * 
 */
void
DR_test();

#endif // !DRICH_EDIT_H

