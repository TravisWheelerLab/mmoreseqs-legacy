/*******************************************************************************
 *  FILE:      hmm_fasta_parser.c
 *  PURPOSE:   Converts a .fasta to HMM_PROFILE.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/* local imports */
#include "objects/structs.h"
#include "objects/hmm_profile.h"
#include "objects/sequence.h"
#include "utility.h"

/* header */
#include "hmm_parser.h"

/* Converts .fasta file to build a HMM_PROFILE object */
HMM_PROFILE* HMM_PROFILE_Fasta_Parse( char*   _filename_,
                                      long    offset )
{
   
}  

/* Converts SEQUENCE object a HMM_PROFILE object */
void SEQUENCE_to_HMM_PROFILE_Convert( SEQUENCE*     seq,
                                      HMM_PROFILE*  prof )
{
   
}  
