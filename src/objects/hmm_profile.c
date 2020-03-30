/*******************************************************************************
 *  @file hmm_profile.c
 *  @brief HMM_PROFILE Object
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

/* local imports */
#include "structs.h"
#include "../utility.h"

/* header */
#include "hmm_profile.h"

/* Constructor */
HMM_PROFILE* HMM_PROFILE_Create()
{
   HMM_PROFILE*   prof = NULL;

   prof = (HMM_PROFILE*) malloc( sizeof(HMM_PROFILE) );
   if (prof == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc HMM_PROFILE.\n");
      exit(EXIT_FAILURE);
   }

   prof->filepath       = NULL;
   prof->name           = NULL;
   prof->acc            = NULL;
   prof->desc           = NULL;
   prof->alph           = NULL;
   
   prof->msv_dist       = NULL;
   prof->viterbi_dist   = NULL;
   prof->forward_dist   = NULL;

   prof->N              = 0;
   prof->Nalloc         = 0;
   prof->alph_leng      = 20;

   prof->bg_model       = NULL;
   prof->hmm_model      = NULL;

   prof->bg_model       = (HMM_BG*) calloc( 1, sizeof(HMM_BG) );
   if (prof->bg_model == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc BG_MODEL in HMM_PROFILE.\n");
      exit(EXIT_FAILURE);
   }

   return prof;
}

/* Destructor */
void HMM_PROFILE_Destroy( HMM_PROFILE* prof )
{
   free(prof->filepath);
   free(prof->name);
   free(prof->acc);
   free(prof->desc);
   free(prof->alph);

   free(prof->msv_dist);
   free(prof->viterbi_dist);
   free(prof->forward_dist);

   free(prof->bg_model);
   free(prof->hmm_model);

   free(prof);
}

/* Set text to given HMM_PROFILE field */
void HMM_PROFILE_Set_TextField( char** prof_field, 
                                char*  text )
{
   *prof_field = malloc( sizeof(char) * ( strlen(text) + 1 ) );
   if ( *prof_field == NULL ) {
      fprintf(stderr, "ERROR: Unable to malloc TEXTFIELD for HMM_PROFILE.\n");
      exit(EXIT_FAILURE);
   }
   strcpy( *prof_field, text );
}

/* Set HMM Model Length and allocate memory for nodes */
void HMM_PROFILE_Set_Model_Length( HMM_PROFILE* prof, 
                                   int          length )
{
   /* realloc memory if allocated length is less than new length */
   if ( prof->Nalloc < length )
   {
      prof->hmm_model = (HMM_NODE*) realloc( prof->hmm_model, (length + 1) * sizeof(HMM_NODE) );
      prof->Nalloc = length;

      if (prof->hmm_model == NULL) {
         fprintf(stderr, "ERROR: Unable to malloc HMM_MODEL for HMM_PROFILE.\n");
         exit(EXIT_FAILURE);
      }
   }

   prof->N = length;
}

/* Set alphabet (DNA or AMINO ACID) for HMM_PROFILE */
void HMM_PROFILE_Set_Alphabet( HMM_PROFILE* prof, 
                               char*        alph_name )
{
   prof->alph = malloc( sizeof(char) * (strlen(alph_name) + 1) );
   if (prof->alph == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc ALPHABET for HMM_PROFILE.\n");
      exit(EXIT_FAILURE);
   }

   int i = 0;
   while (alph_name[i] != '\0') {
      alph_name[i] = tolower(alph_name[i]);
      i++;
   }

   if ( strcmp(alph_name, "amino") == 0 ) {
      prof->alph_leng = NUM_AMINO;
   }
   else if ( strcmp(alph_name, "dna") == 0 ) {
      prof->alph_leng = NUM_DNA;
   }
   else {
      fprintf(stderr, "ERROR: Invalid alphabet type: %s\n", alph_name );
      exit(EXIT_FAILURE);
   }
}

/* Set Distribution Parameters for HMM_PROFILE */
void HMM_PROFILE_Set_Distribution_Params( HMM_PROFILE* prof, 
                                          float        param1, 
                                          float        param2, 
                                          char*        dist_name )
{
   float*   parPtr1 = NULL;
   float*   parPtr2 = NULL;

   if ( strcmp( dist_name, "MSV" ) == 0 ) {
      parPtr1 = &prof->msv_dist->param1;
      parPtr2 = &prof->msv_dist->param2;
   }
   else if ( strcmp( dist_name, "VITERBI" ) == 0 ) {
      parPtr1 = &prof->viterbi_dist->param1;
      parPtr2 = &prof->viterbi_dist->param2;
   }
   else if ( strcmp( dist_name, "FORWARD" ) == 0 ) {
      parPtr1 = &prof->forward_dist->param1;
      parPtr2 = &prof->forward_dist->param2;
   }
   else {
      fprintf( stderr, "Invalid distribution type: %s\n", dist_name );
      exit(EXIT_FAILURE);
   }

   parPtr1 = &param1;
   parPtr2 = &param2;
}

/* Output HMM_PROFILE to FILE POINTER */
void HMM_PROFILE_Dump( HMM_PROFILE* prof, 
                       FILE*        fp )
{
   /* test for bad file pointer */
   if (fp == NULL) {
      const char* obj_name = "HMM_PROFILE";
      fprintf(stderr, "ERROR: Bad FILE POINTER for printing %s.\n", obj_name);
      exit(EXIT_FAILURE);
      return;
   }

   int i,j;
   float best_val, new_val;
   char best_amino;
   HMM_NODE curr_node;

   fprintf(fp, "\n");
   fprintf(fp, "===== HMM PROFILE ====================================\n");
   fprintf(fp, "\t         NAME:\t%s\n", prof->name);
   fprintf(fp, "\t       LENGTH:\t%d\n", prof->N);
   char *seq = (char *)malloc( sizeof(char) * prof->N );

   /* print consensus sequence */
   fprintf(fp, "\tCONSENSUS SEQ:\t");
   for (int i = 1; i < prof->N+1; i++)
   {
      best_val = -INF;
      curr_node = prof->hmm_model[i];
      for (int j = 0; j < NUM_AMINO; j++)
      {
         new_val = curr_node.match[j];
         if (best_val < new_val)
         {
            best_val = new_val;
            best_amino = AA[j];
         }
      }
      fprintf(fp, "%c", best_amino);
   }
   fprintf(fp, "\n");

   fprintf(fp, "\t        COMPO:\n");
   fprintf(fp, "FREQ:\t");
   for (int j = 0; j < NUM_AMINO; j++)
   {
      fprintf(fp, "\t%.4f", prof->bg_model->freq[j]);
   }
   fprintf(fp, "\nCOMPO:\t");
   for (int j = 0; j < NUM_AMINO; j++)
   {
      fprintf(fp, "\t%.4f", prof->bg_model->compo[j]);
   }
   fprintf(fp, "\nINSERT:\t");
   for (int j = 0; j < NUM_AMINO; j++)
   {
      fprintf(fp, "\t%.4f", prof->bg_model->insert[j]);
   }
   fprintf(fp, "\nTRANS:\t");
   for (int j = 0; j < NUM_TRANS_STATES; j++)
   {
      fprintf(fp, "\t%.4f", prof->bg_model->trans[j]);
   }
   fprintf(fp, "\n\n");

   /* NORMAL STATE PROBS */
   fprintf(fp, "=== TRANSITION PROBS ===\n");
   for (int i = 0; i < prof->N+1; i++)
   {
      fprintf(fp, "%d", i);
      for (int j = 0; j < NUM_TRANS_STATES; j++)
      {
         fprintf(fp, "\t%.4f", prof->hmm_model[i].trans[j]);
      }
      fprintf(fp, "\n");
   }
   fprintf(fp, "=== MATCH EMISSION PROBS ===\n");
   for (int i = 0; i < prof->N+1; i++)
   {
      fprintf(fp, "%d", i);
      for (int j = 0; j < NUM_AMINO; j++)
      {
         fprintf(fp, "\t%.4f", prof->hmm_model[i].match[j]);
      }
      fprintf(fp, "\n");
   }
   fprintf(fp, "=== INSERT EMISSION PROBS ===\n");
   for (int i = 0; i < prof->N+1; i++)
   {
      fprintf(fp, "%d", i);
      for (int j = 0; j < NUM_AMINO; j++)
      {
         fprintf(fp, "\t%.4f", prof->hmm_model[i].insert[j]);
      }
      fprintf(fp, "\n");
   }

   fprintf(fp, "BACKGROUND:\n");
   fprintf(fp, "LOG:\t\t");
   for (int i = 0; i < NUM_AMINO; i++)
   {
      fprintf(fp, "%.4f\t", BG_MODEL_log[i]);
   }
   fprintf(fp, "\n");
   fprintf(fp, "ACTUAL:\t\t");
   for (int i = 0; i < NUM_AMINO; i++)
   {
      fprintf(fp, "%.4f\t", BG_MODEL[i]);
   }
   fprintf(fp, "\n\n");

   fprintf(fp, "SPECIAL:\n");
   fprintf(fp, "\tE:\t%.4f\t%.4f\n", prof->bg_model->spec[SP_E][SP_LOOP], prof->bg_model->spec[SP_E][SP_MOVE]);
   fprintf(fp, "\tN:\t%.4f\t%.4f\n", prof->bg_model->spec[SP_N][SP_LOOP], prof->bg_model->spec[SP_N][SP_MOVE]);
   fprintf(fp, "\tC:\t%.4f\t%.4f\n", prof->bg_model->spec[SP_C][SP_LOOP], prof->bg_model->spec[SP_C][SP_MOVE]);
   fprintf(fp, "\tJ:\t%.4f\t%.4f\n", prof->bg_model->spec[SP_J][SP_LOOP], prof->bg_model->spec[SP_J][SP_MOVE]);

   fprintf(fp, "//");
}

