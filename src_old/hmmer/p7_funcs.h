/*******************************************************************************
 *  FILE:      p7_funcs.h
 *  PURPOSE:   Imports HMMER and Easel library functions.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _P7_FUNCS_H_
#define _P7_FUNCS_H_

P7_BG *
p7_bg_Create(const ESL_ALPHABET *abc);

P7_BUILDER *
p7_builder_Create(const ESL_GETOPTS *go, const ESL_ALPHABET *abc);

void
p7_bg_Destroy(P7_BG *bg);

int
p7_SingleBuilder( P7_BUILDER*    bld, 
                  ESL_SQ*        sq, 
                  P7_BG*         bg, 
                  P7_HMM**       opt_hmm,
                  P7_TRACE**     opt_tr, 
                  P7_PROFILE**   opt_gm, 
                  P7_OPROFILE**  opt_om );

int
p7_Seqmodel( const ESL_ALPHABET* abc, 
                   ESL_DSQ*      dsq, 
                   int           M, 
                   char*         name,
                   ESL_DMATRIX*  Q, 
                   float*        f, 
                   double        popen, 
                   double        pextend,
                   P7_HMM**      ret_hmm);

int
p7_hmm_SetComposition(P7_HMM *hmm);

int
p7_hmm_SetConsensus(P7_HMM*   hmm, 
                    ESL_SQ*   sq);

int
p7_hmm_CalculateOccupancy(const P7_HMM*   hmm, 
                                float*    mocc, 
                                float*    iocc);

#endif /* _P7_FUNCS_H_ */

