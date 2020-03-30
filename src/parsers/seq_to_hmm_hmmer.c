/*******************************************************************************
 *  FILE:      seq_to_model.c
 *  PURPOSE:   Converts a SEQUENCE to HMM_PROFILE.
 *             Uses the Easel and HMMER libraries.
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
#include "objects/score_matrix.h"
#include "utility.h"

/* hmmer imports */
// #include "hmmer.h"

/* self header */
#include "seq_to_model.h"

/* converts SEQUENCE to HMM_PROFILE model */
void SEQUENCE_to_HMM_PROFILE( SEQUENCE*      seq, 
                              HMM_PROFILE*   prof )
{

}
