/*******************************************************************************
 *  FILE:      parser.h
 *  PURPOSE:   Parses .hmm files into HMM_PROFILE object
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _HMM_PARSER_H
#define _HMM_PARSER_H

/* === INCLUDES === */
// #include "objects/structs.c"
// #include "objects/sequence.c"
// #include "objects/hmm_profile.c"

/* === FUNCTIONS === */

/* parse .hmm file and build HMM_PROFILE object */
void HMM_PROFILE_Parse( HMM_PROFILE*   prof,
                        char*          _filename_,
                        long           offset );

void HMM_PROFILE_Parse_REAL(  HMM_PROFILE*   prof,
                              char*          _filename_,
                              long           offset );

/* .hmm stores numbers in log space, but we need reals */
void HMM_PROFILE_Convert_NegLog_To_Real( HMM_PROFILE* prof );

/* NOTE: modeled after HMMER p7_ProfileConfig() and other p7 functions */
/* Configures HMM_PROFILE to account for background model */
void HMM_PROFILE_Config( HMM_PROFILE* prof, 
                         int          mode );
/* Calculates the Occupancy for the HMM_PROFILE */
void HMM_PROFILE_CalcOccupancy( HMM_PROFILE* prof, 
                                float*       mocc,
                                float*       iocc );
/* Reconfigure the Length of the HMM_PROFILE */
void HMM_PROFILE_ReconfigLength( HMM_PROFILE*  prof, 
                                 int           L );

#endif /* _HMM_PARSER_H */
