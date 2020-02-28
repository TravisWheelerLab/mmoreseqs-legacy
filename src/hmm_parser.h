/*******************************************************************************
 *  @file parser.h
 *  @brief Parses .hmm and .fasta files
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

#ifndef _HMM_PARSER_H
#define _HMM_PARSER_H

/* === INCLUDES === */
// #include "objects/structs.c"
// #include "objects/sequence.c"
// #include "objects/hmm_profile.c"

/* === FUNCTIONS === */
/* parse .hmm file and build HMM_PROFILE object */
HMM_PROFILE* HMM_PROFILE_Parse(char* _filename_);

/* configures HMM_PROFILE to account for background model */
/* modeled after HMMER p7_ProfileConfig() */
void HMM_PROFILE_Config(HMM_PROFILE* prof, 
                        int          mode);
/* Calculates the Occupancy for the HMM_PROFILE */
void HMM_PROFILE_CalcOccupancy(HMM_PROFILE* prof, 
                               float*       occ);
/* Reconfigure the Length of the HMM_PROFILE */
void HMM_PROFILE_ReconfigLength(HMM_PROFILE*  prof, 
                                int           L);

/* parse .fasta file and build SEQUENCE object */
SEQUENCE* SEQUENCE_Parse(char*  _filename_);

#endif /* _HMM_PARSER_H */
