/*******************************************************************************
 *  FILE:      seq_to_model.c
 *  PURPOSE:   Converts a SEQUENCE to HMM_PROFILE.
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

/* self header */
#include "seq_to_model.h"

/* converts SEQUENCE to HMM_PROFILE model */
void SEQ_to_PROF( SEQUENCE*      seq, 
                  HMM_PROFILE*   prof )
{
   char*          fname          = NULL;     /* filename of sequence */
   int            N              = 0;        /* length of sequence */

   double         slambda        = 0.;       /* */
   double         popen          = 0.;       /* gap open probability */
   double         pextend        = 0.;       /* gap extention probability */

   char*          submat_fname   = NULL;     /* substitution matrix filename */
   SCORE_MATRIX*  submat         = NULL;     /* substitution matrix */
   HMM_BG*        hmm_bg         = NULL;     /* background null model */

   /* gap open/extend (defaults from HMMER) */
   popen    = 0.02;
   pextend  = 0.4;

   /* substitution matrix used to assess seq */
   submat_fname   = "../data/submat/blosum62.submat";
   submat         = SCORE_MATRIX_Load(submat_fname);

   /* background model using hardcoded frequencies */
   hmm_bg         = HMM_BG_Create();

   /* TODO: complete joint probability conversion */
   /* calculate joint probability */
   // SCORE_MATRIX_Prob_Given_BG( submat, &slambda, & )

   /* extract info */
   N        = seq->N;
   fname    = seq->filename;

   /* resize profile if necessary */
   HMM_PROFILE_Set_Model_Length( prof, seq->N );

   for ( int i = 0; i < N; i++ )
   {
      /* match emission */
      for ( int j = 0; j < NUM_AMINO; j++ ) {

      }

      /* insertion emmission */
      for ( int j = 0; j < NUM_AMINO; j++ ) {

      }

      /* transition scores */
      TSC_HMM( prof, i, M2M ) = 1.0 - popen;
      TSC_HMM( prof, i, M2I ) = popen;
      TSC_HMM( prof, i, M2D ) = popen;
      TSC_HMM( prof, i, I2M ) = 1.0 - pextend;
      TSC_HMM( prof, i, I2I ) = pextend;
      TSC_HMM( prof, i, D2M ) = 1.0 - pextend;
      TSC_HMM( prof, i, D2D ) = pextend;

   }

   /* Final node */
   TSC_HMM( prof, N, M2M ) = 1.0 - popen;
   TSC_HMM( prof, N, M2D ) = 0.0;
   TSC_HMM( prof, N, D2M ) = 1.0;
   TSC_HMM( prof, N, D2D ) = 0.0;
}

HMM_BG* HMM_BG_Create()
{
   HMM_BG* bg = (HMM_BG*) malloc( sizeof(HMM_BG) );
   if (bg == NULL) {
      fprintf( stderr, "ERROR: Unable to malloc for HMM_BG.\n" );
      exit(EXIT_FAILURE);
   }

   /* set background frequencies from BG_MODEL (imported from HMMER) */
   for ( int i = 0; i < NUM_AMINO; i++ ) {
      bg->freq[i] = BG_MODEL[i];
   }

   return bg;
}



