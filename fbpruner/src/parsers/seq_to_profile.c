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
#include "../objects/structs.h"
#include "../utilities/_utilities.h"
#include "../objects/_objects.h"

/* self header */
#include "_parsers.h"

/* converts single SEQUENCE to HMM_PROFILE model */
void SEQUENCE_to_HMM_PROFILE( SEQUENCE*      seq, 
                              HMM_PROFILE*   prof )
{
   printf("# converting sequence to profile...\n");
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
   char           a              = 0;        /* profile emission */
   char           b              = 0;        /* sequence emission */

   double         popen          = 0.;       /* gap open probability */
   double         pextend        = 0.;       /* gap extention probability */

   char*          bld_path       = NULL;     /* path to data folder */
   char*          bld_fname      = NULL;     /* path to singlebuilder matrix */
   /* NOTE: this is global so it only needs to be initialized once */
   // SCORE_MATRIX*  bld            = NULL;

   /* gap open/extend (defaults from HMMER) */
   popen    = 0.02;
   pextend  = 0.4;

   /* SingleBuilder probability matrix (created by HMMER using p7_SingleBuilder, dependent on only bg frequencies) */
   /* NOTE: currently stored externally, but could be integrated in worker object.  Also, never freed. */
   /* TODO: need to set PREFIX in makefile */
   if ( bld == NULL ) {
      bld_path       = strdup( MACRO_STR(PREFIX) );
      bld_fname      = strdup( "data/submat/singlebuilder.submat" );
      bld            = SCORE_MATRIX_Load( bld_fname );
      ERROR_free( bld_path );
      ERROR_free( bld_fname );
   }

   /* resize profile if necessary */
   N  = seq->N;
   HMM_PROFILE_Set_Model_Length( prof, seq->N );
   HMM_PROFILE_Set_TextField( &prof->name, seq->name );

   /* special match probabilities for initial node */
   MSC_X( prof, 0, 0 ) = 1.;
   for ( int j = 1; j < NUM_AMINO; j++ )
      MSC_X( prof, 0, j ) = 0.;

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
            MSC_X( prof, i, j ) = SCORE_MATRIX_Get_Score( bld, a, b );
         }
      }

      /* insertion emmission (uses hardcoded background frequencies) */
      for ( int j = 0; j < NUM_AMINO; j++ ) {
         ISC_X( prof, i, j ) = BG_MODEL[j];
      }

      /* transition scores */
      TSC_X( prof, i, M2M ) = 1.0 - 2 * popen;
      TSC_X( prof, i, M2I ) = popen;
      TSC_X( prof, i, M2D ) = popen;
      TSC_X( prof, i, I2M ) = 1.0 - pextend;
      TSC_X( prof, i, I2I ) = pextend;
      TSC_X( prof, i, D2M ) = 1.0 - pextend;
      TSC_X( prof, i, D2D ) = pextend;
   }

   /* final node transitions */
   TSC_X( prof, N, M2M ) = 1.0 - popen;
   TSC_X( prof, N, M2D ) = 0.0;
   TSC_X( prof, N, D2M ) = 1.0;
   TSC_X( prof, N, D2D ) = 0.0;

   /* update bg data to reflect insert and transition (same as all but first and last node) */
   for ( int j = 0; j < NUM_AMINO; j++ ) {
      prof->bg_model->insert[j] = prof->hmm_model[1].insert[j];
   }
   for ( int j = 0; j < NUM_TRANS_STATES; j++ )  {
      prof->bg_model->trans[j] = prof->hmm_model[1].trans[j];
   }

   /* free data */
   // SCORE_MATRIX_Destroy( bld );
}

/* set background composition of sequence */
void HMM_PROFILE_Set_Composition( HMM_PROFILE* prof )
{
   float*   mocc = NULL;   /* match occupancy */
   float*   iocc = NULL;   /* insert occupancy */
   int      k;

   /* */
   mocc = (float*) ERROR_malloc( sizeof(float) * (prof->N + 1) );
   iocc = (float*) ERROR_malloc( sizeof(float) * (prof->N + 1) );

   /* calculate occupancy (counting) to update composition */
   HMM_PROFILE_CalcOccupancy( prof, mocc, iocc );
   for (int i = 0; i < NUM_AMINO; i++) {
      prof->bg_model->compo[i] = 0.0;
   }
   for (int i = 0; i < NUM_AMINO; i++) {
      prof->bg_model->compo[i] += prof->hmm_model[0].insert[i] * iocc[0];
   }
   for (int j = 1; j <= prof->N; j++) 
   {
      for (int i = 0; i < NUM_AMINO; i++)
         prof->bg_model->compo[i] += prof->hmm_model[j].match[i] * mocc[j];
      for (int i = 0; i < NUM_AMINO; i++)
         prof->bg_model->compo[i] += prof->hmm_model[j].insert[i] * iocc[j];
   }

   /* normalize the composition */
   int sum = 0;
   for (int i = 0; i < NUM_AMINO; i++) {
      sum += prof->bg_model->compo[i];
   }
   if (sum != 0.0) {
      for (int i = 0; i < NUM_AMINO; i++) {
         prof->bg_model->compo[i] /= sum;
      }
   } else {
      for (int i = 0; i < NUM_AMINO; i++) {
         prof->bg_model->compo[i] = 1 / (float) NUM_AMINO;
      }
   }

   ERROR_free(mocc);
   ERROR_free(iocc);
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
   // ESL_RANDOMNESS*   r        = NULL;
   // int               EmL      = 200;
   // int               EmN      = 200;
   // int               EvL      = 200;
   // int               EvN      = 200;
   // int               EfL      = 100;
   // int               EfN      = 200;
   // double            Eft      = 0.04;
   // double            lambda, mmu, vmu, tau;

   // /* initialize random number generator */
   // ((r = esl_randomness_CreateFast(42)) == NULL) ESL_XFAIL(eslEMEM, errbuf, "failed to create RNG");
   // /* lambda: MSV, VITERBI, FORWARD */
   // HMM_PROFILE_CalibrateLambda( prof, &lambda );
   // /* mu: MSV */
   // // HMM_PROFILE_CalibrateMSVMu( prof, r, EmL, EmN, lambda, &mmu );
   // /* mu: Viterbi */
   // // HMM_PROFILE_CalibrateViterbiMu( prof, r, EmL, EmN, lambda, &mmu );
   // /* tau */
   // HMM_PROFILE_CalibrateTau( prof, r, EfL, EfN, lambda, Eft, &tau );

   // esl_randomness_Destroy(r);
}

// /* modeled after: p7_Lambda() */
// void HMM_PROFILE_CalibrateLambda(   HMM_PROFILE*   prof, 
//                                     int*           ret_lambda )
// {
//    double H = HMM_PROFILE_MeanMatchRelativeEntropy( prof );

//    *ret_lambda = eslCONST_LOG2 + 1.44 / ((double) hmm->M * H);
// }

// /* modeled after: p7_MeanMatchRelativeEntropy() */
// double HMM_PROFILE_MeanMatchRelativeEntropy( HMM_PROFILE* prof )
// {
//    int         k;
//    double      KL = 0.0;
//    HMM_COMPO*  bg = prof->bg_model;

//    for (k = 1; k <= prof->N; k++) {
//       KL += esl_vec_FRelEntropy( prof->hmm_model[k].match, bg->freq, NUM_AMINO);
//    }
//    KL /= (double) prof->N;
//    return KL;
// }

// /* modeled after: esl_vec_FRelEntropy() */
// float HMM_PROFILE_RelativeEntropy(const float *p, const float *q, int n)
// {
//   int    i;
//   float  kl;

//   kl = 0.;
//   for(i = 0; i < n; i++)
//     if (p[i] > 0.) {
//       if (q[i] == 0.) return eslINFINITY;
//       else            kl += p[i] * log2(p[i]/q[i]);
//     }
//   return kl;
// }



