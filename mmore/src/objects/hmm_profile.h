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
HMM_PROFILE* 
HMM_PROFILE_Create();

/* Destructor */
HMM_PROFILE* 
HMM_PROFILE_Destroy( HMM_PROFILE*   prof );

/* reuse profile by setting length of length of profile to zero */
int 
HMM_PROFILE_Reuse( HMM_PROFILE*  prof );

/* Create backround hmm from hardcoded background frequencies */
HMM_COMPO* 
HMM_COMPO_Create();

/* Set Textfield to HMM_PROFILE field */
void 
HMM_PROFILE_SetTextField( char**   prof_field, 
                           char*    text );

/* Set HMM Model Length and allocate memory for nodes */
void 
HMM_PROFILE_SetModel_Length( HMM_PROFILE*    prof, 
                              int             length );

/* Set alphabet (DNA or AMINO ACID) for HMM_PROFILE */
void 
HMM_PROFILE_SetAlphabet(  HMM_PROFILE*     prof, 
                           char*            alph_name );

/* Determine the consensus sequence (highest likelihood output) using HMM_PROFILE */
void 
HMM_PROFILE_SetConsensus( HMM_PROFILE*    prof );

/**  FUNCTION:  HMM_PROFILE_GetConsensus()
 *   SYNOPSIS:  Return <consensus> sequence as string.
 *              Will build it now if necessary.
 */
STR 
HMM_PROFILE_GetConsensus( HMM_PROFILE*    prof );

/* Set Distribution Parameters for HMM_PROFILE */
void 
HMM_PROFILE_SetDistribution_Params(   HMM_PROFILE*   prof, 
                                       float          param1, 
                                       float          param2, 
                                       char*          dist_name );

/* modeled after p7_bg_SetLength()  */
void HMM_PROFILE_BG_SetLength(   HMM_PROFILE*   prof,
                                 int            length );

/* Constrain model to range [tbeg, tend] */
void 
HMM_PROFILE_SetSubmodel( HMM_PROFILE* prof, int tbeg, int tend );

/* Unconstrain model to cover entire */
void 
HMM_PROFILE_UnsetSubmodel( HMM_PROFILE* prof );

/* Output HMM_PROFILE to FILE POINTER */
void 
HMM_PROFILE_Dump( HMM_PROFILE*    prof, 
                  FILE*           fp );

#endif /* _HMM_PROFILE_H */
