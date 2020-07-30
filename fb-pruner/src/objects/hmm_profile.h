/*******************************************************************************
 *  FILE:      hmm_profile.c
 *  PURPOSE:   HMM_PROFILE Object.
 *
 *  AUTHOR:    Dave Rich
 *  BUGS:  
 *******************************************************************************/

#ifndef _HMM_PROFILE_H
#define _HMM_PROFILE_H

/* Constructor */
HMM_PROFILE* HMM_PROFILE_Create();

/* Destructor */
void* HMM_PROFILE_Destroy( HMM_PROFILE* prof );

/* reuse profile by setting length of length of profile to zero */
void HMM_PROFILE_Reuse( HMM_PROFILE* prof );

/* Create backround hmm from hardcoded background frequencies */
HMM_COMPO* HMM_COMPO_Create();

/* Set Textfield to HMM_PROFILE field */
void HMM_PROFILE_Set_TextField( char** prof_field, 
                                char*  text );

/* Set HMM Model Length and allocate memory for nodes */
void HMM_PROFILE_Set_Model_Length( HMM_PROFILE* prof, 
                                   int          length );

/* Set alphabet (DNA or AMINO ACID) for HMM_PROFILE */
void HMM_PROFILE_Set_Alphabet( HMM_PROFILE* prof, 
                               char*        alph_name );

/* Determine the consensus sequence (highest likelihood output) using HMM_PROFILE */
void HMM_PROFILE_Set_Consensus( HMM_PROFILE* prof );

/* Set Distribution Parameters for HMM_PROFILE */
void HMM_PROFILE_Set_Distribution_Params( HMM_PROFILE* prof, 
                                          float        param1, 
                                          float        param2, 
                                          char*        dist_name );

/* Output HMM_PROFILE to FILE POINTER */
void HMM_PROFILE_Dump( HMM_PROFILE* prof, 
                       FILE*        fp );

#endif /* _HMM_PROFILE_H */
