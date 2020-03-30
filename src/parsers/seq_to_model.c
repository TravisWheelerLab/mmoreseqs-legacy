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
#include "utility.h"
/* self header */
#include "seq_to_model.h"

/* converts sequence to HMM profile model */
void SEQ_to_PROF( SEQUENCE*      seq, 
                  HMM_PROFILE*   prof )
{
   int      N;
   char*    filename;
   double   popen;
   double   pextend;

   /* extract info */
   N        = seq->N;
   filename = seq->filename;

   /* resize profile if necessary */
   HMM_PROFILE_Set_Model_Length( prof, seq->N );

   for (int i = 0; i < N; i++)
   {
      for (int j = 0; j < NUM_AMINO; j++) {

      }
   }
}