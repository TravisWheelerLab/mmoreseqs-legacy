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
#include "structs.h"
#include "utilities.h"
#include "objects.h"

/* self header */
#include "parsers.h"

/* converts single SEQUENCE to HMM_PROFILE model */
void SEQUENCE_to_HMM_PROFILE( SEQUENCE*      seq, 
                              HMM_PROFILE*   prof )
{
   printf("converting sequence to profile...\n");
   HMM_PROFILE_From_Seq( prof, seq );
   HMM_PROFILE_Set_Composition( prof );
   HMM_PROFILE_Calibrate( prof );
}

/* create HMM from SEQUENCE */
void HMM_PROFILE_From_Seq( HMM_PROFILE*   prof,
                           SEQUENCE*      seq )
{
   FILE*          fp             = NULL;     /* debugger filepointer */

   int            N              = 0;        /* length of sequence */
   char           a              = 0;        /* target emission */
   char           b              = 0;        /* query emission */

   double         popen          = 0.;       /* gap open probability */
   double         pextend        = 0.;       /* gap extention probability */

   char*          bld_fname      = NULL;     /* path to singlebuilder matrix */
   SCORE_MATRIX*  bld            = NULL;     /* SingleBuilder probability matrix (created by HMMER using p7_SingleBuilder, dependent on only bg frequencies) */

   /* gap open/extend (defaults from HMMER) */
   popen    = 0.02;
   pextend  = 0.4;

   /* singlebuilder matrix */
   /* NOTE: could let this persist when creating multiple sequence-to-profiles */
   bld_fname      = "data/submat/singlebuilder.submat";
   bld            = SCORE_MATRIX_Load( bld_fname );

   /* resize profile if necessary */
   N  = seq->N;
   HMM_PROFILE_Set_Model_Length( prof, seq->N );
   HMM_PROFILE_Set_TextField( &prof->name, seq->name );

   /* special match probabilities for initial node */
   MSC_HMM( prof, 0, 0 ) = 1.;
   for ( int j = 1; j < NUM_AMINO; j++ )
      MSC_HMM( prof, 0, j ) = 0.;

   /* for each node in sequence */
   for ( int i = 0; i <= N; i++ )
   {
      if (i > 0) 
      {
         a = seq->seq[i - 1];

         /* match emission (uses singlebuilder ) */
         for ( int j = 0; j < NUM_AMINO; j++ ) 
         {
            b = AA[j];
            MSC_HMM( prof, i, j ) = SCORE_MATRIX_Get_Score( bld, a, b );
         }
      }

      /* insertion emmission (uses hardcoded background frequencies) */
      for ( int j = 0; j < NUM_AMINO; j++ ) 
      {
         ISC_HMM( prof, i, j ) = BG_MODEL[j];
      }

      /* transition scores */
      TSC_HMM( prof, i, M2M ) = 1.0 - 2 * popen;
      TSC_HMM( prof, i, M2I ) = popen;
      TSC_HMM( prof, i, M2D ) = popen;
      TSC_HMM( prof, i, I2M ) = 1.0 - pextend;
      TSC_HMM( prof, i, I2I ) = pextend;
      TSC_HMM( prof, i, D2M ) = 1.0 - pextend;
      TSC_HMM( prof, i, D2D ) = pextend;
   }

   /* final node transitions */
   TSC_HMM( prof, N, M2M ) = 1.0 - popen;
   TSC_HMM( prof, N, M2D ) = 0.0;
   TSC_HMM( prof, N, D2M ) = 1.0;
   TSC_HMM( prof, N, D2D ) = 0.0;

   /* update bg data to reflect insert and transition (same as all but first and last node) */
   for ( int j = 0; j < NUM_AMINO; j++ )
      prof->bg_model->insert[j] = prof->hmm_model[1].insert[j];
   for ( int j = 0; j < NUM_TRANS_STATES; j++ ) 
      prof->bg_model->trans[j] = prof->hmm_model[1].trans[j];

   /* free data */
   SCORE_MATRIX_Destroy( bld );
}

/* set background composition of sequence */
void HMM_PROFILE_Set_Composition( HMM_PROFILE* prof )
{
   float*   mocc = NULL;   /* match occupancy */
   float*   iocc = NULL;   /* insert occupancy */
   int      k;

   /* */
   mocc = (float*) malloc( sizeof(float) * (prof->N + 1) );
   if ( mocc == NULL ) {
      fprintf(stderr, "ERROR: Unable to malloc.\n");
      exit(EXIT_FAILURE);
   }
   iocc = (float*) malloc( sizeof(float) * (prof->N + 1) );
   if ( iocc == NULL ) {
      fprintf(stderr, "ERROR: Unable to malloc.\n");
      exit(EXIT_FAILURE);
   }

   /* calculate occupancy (counting) to update composition */
   HMM_PROFILE_CalcOccupancy( prof, mocc, iocc );
   for (int i = 0; i < NUM_AMINO; i++) 
      prof->bg_model->compo[i] = 0.0;
   for (int i = 0; i < NUM_AMINO; i++)
      prof->bg_model->compo[i] += prof->hmm_model[0].insert[i] * iocc[0];
   for (int j = 1; j <= prof->N; j++) 
   {
      for (int i = 0; i < NUM_AMINO; i++)
         prof->bg_model->compo[i] += prof->hmm_model[j].match[i] * mocc[j];
      for (int i = 0; i < NUM_AMINO; i++)
         prof->bg_model->compo[i] += prof->hmm_model[j].insert[i] * iocc[j];
   }

   /* normalize the composition */
   int sum = 0;
   for (int i = 0; i < NUM_AMINO; i++)
      sum += prof->bg_model->compo[i];
   if (sum != 0.0) {
      for (int i = 0; i < NUM_AMINO; i++)
         prof->bg_model->compo[i] /= sum;
   } else {
      for (int i = 0; i < NUM_AMINO; i++)
         prof->bg_model->compo[i] = 1 / (float) NUM_AMINO;
   }

   free(mocc);
   free(iocc);
}

/* TODO */
// void HMM_BG_Set_Length( HMM_BG*  bg, 
//                         int      L )
// {
//    bg->p1 = (float) L / (float) (L+1);

//    bg->fhmm->t[0][0] = bg->p1;
//    bg->fhmm->t[0][1] = 1.0f - bg->p1;
// }

/* TODO */
/* run simulation to calibrate e-value parameters of model */
void HMM_PROFILE_Calibrate( HMM_PROFILE* prof )
{

}


