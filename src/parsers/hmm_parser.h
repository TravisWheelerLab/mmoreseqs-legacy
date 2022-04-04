/*******************************************************************************
 *  - FILE:   parser.h
 *  - DESC:    Parses .hmm files into HMM_PROFILE object
 *******************************************************************************/

#ifndef _HMM_PARSER_H
#define _HMM_PARSER_H

/* === FUNCTIONS === */

/* parse .hmm file and build HMM_PROFILE object */
void HMM_PROFILE_Parse(HMM_PROFILE* prof, char* filename, long offset);

/* computes the expected value for the match state of current node */
float HMM_NODE_Expected_Value(HMM_NODE* node);

/* .hmm stores numbers in log space, but we need reals */
void HMM_PROFILE_Convert_NegLog_To_Real(HMM_PROFILE* prof);

/* NOTE: modeled after HMMER p7_ProfileConfig() and other p7 functions */
/* Configures HMM_PROFILE to account for background model */
void HMM_PROFILE_Config(HMM_PROFILE* prof, int mode);

/* Calculates the Occupancy for the HMM_PROFILE */
void HMM_PROFILE_CalcOccupancy(HMM_PROFILE* prof, float* mocc, float* iocc);

/* Reconfigure the Length of the HMM_PROFILE */
void HMM_PROFILE_ReconfigLength(HMM_PROFILE* prof, int L);

/* Configure the Length of the HMM_PROFILE based on the length of the sequence
 */
void HMM_PROFILE_ReconfigUnihit(HMM_PROFILE* prof, int L);

/* Configure the Length of the HMM_PROFILE based on the length of the sequence
 */
void HMM_PROFILE_ReconfigMultihit(HMM_PROFILE* prof, int L);

#endif /* _HMM_PARSER_H */
