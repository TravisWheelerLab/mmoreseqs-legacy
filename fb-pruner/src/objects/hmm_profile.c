/*******************************************************************************
 *  FILE:      hmm_profile.c
 *  PURPOSE:   HMM_PROFILE Object.
 *
 *  AUTHOR:    Dave Rich
 *  BUGS:  
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
#include "utilities.h"
#include "objects.h"

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
   prof->consensus      = NULL;
   
   prof->msv_dist       = (DIST_PARAM){0.0, 0.0};
   prof->viterbi_dist   = (DIST_PARAM){0.0, 0.0};
   prof->forward_dist   = (DIST_PARAM){0.0, 0.0};

   prof->N              = 0;
   prof->Nalloc         = 0;
   prof->alph_leng      = 20;

   prof->bg_model       = NULL;
   prof->hmm_model      = NULL;

   prof->bg_model = HMM_COMPO_Create();

   return prof;
}

/* Destructor */
void* HMM_PROFILE_Destroy( HMM_PROFILE* prof )
{
   if (prof == NULL) return prof;
   
   ERROR_free(prof->filepath);
   ERROR_free(prof->name);
   ERROR_free(prof->acc);
   ERROR_free(prof->desc);
   ERROR_free(prof->alph);

   ERROR_free(prof->bg_model);
   ERROR_free(prof->hmm_model);

   ERROR_free(prof);
   prof = NULL;
   return prof;
}

/* reuse profile by setting length of length of profile to zero */
void HMM_PROFILE_Reuse( HMM_PROFILE* prof )
{
   prof->N = 0;
}

/* create backround hmm composition from hardcoded background frequencies */
HMM_COMPO* HMM_COMPO_Create()
{
   HMM_COMPO* bg = (HMM_COMPO*) malloc( sizeof(HMM_COMPO) );
   if (bg == NULL) {
      fprintf( stderr, "ERROR: Unable to malloc for HMM_COMPO.\n" );
      exit(EXIT_FAILURE);
   }

   /* set background frequencies from BG_MODEL (imported from HMMER) */
   for ( int i = 0; i < NUM_AMINO; i++ ) {
      bg->freq[i] = BG_MODEL[i];
   }

   return bg;
}

/* Set text to given HMM_PROFILE field */
void HMM_PROFILE_Set_TextField( char** prof_field, 
                                char*  text )
{
   *prof_field = realloc( *prof_field, sizeof(char) * ( strlen(text) + 1 ) );
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
      prof->hmm_model   = (HMM_NODE*) realloc( prof->hmm_model, (length + 1) * sizeof(HMM_NODE) );
      prof->Nalloc      = length;

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

/* Determine the consensus sequence (highest likelihood output) using HMM_PROFILE */
void HMM_PROFILE_Set_Consensus( HMM_PROFILE* prof )
{
   float       best_val; 
   float       new_val;
   char        best_amino;
   HMM_NODE    curr_node;

   /* TODO: update to allocate in Create() / change to VECTOR_CHAR */
   /* clear pre-existing consensus */
   ERROR_free(prof->consensus);
   /* allocate new consensus */
   prof->consensus = (char*) malloc( sizeof(char) * (prof->N + 1) );

   /* find consensus sequence */
   for (int i = 1; i <= prof->N; i++)
   {
      best_val = -INF;
      curr_node = prof->hmm_model[i];
      /* for each residue in the consensus sequence, choose the most likely to be emitted */
      for (int j = 0; j < NUM_AMINO; j++)
      {
         new_val = curr_node.match[j];
         if (best_val < new_val)
         {
            best_val = new_val;
            best_amino = AA[j];
         }
      }
      prof->consensus[i-1] = best_amino;
   }
   prof->consensus[prof->N] = '\0';
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
      prof->msv_dist.param1 = param1;
      prof->msv_dist.param2 = param2;
   }
   else if ( strcmp( dist_name, "VITERBI" ) == 0 ) {
      prof->viterbi_dist.param1 = param1;
      prof->viterbi_dist.param2 = param2;
   }
   else if ( strcmp( dist_name, "FORWARD" ) == 0 ) {
      prof->forward_dist.param1 = param1;
      prof->forward_dist.param2 = param2;
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

   int i, j;
   if ( prof->consensus == NULL ) {
      HMM_PROFILE_Set_Consensus( prof );
   }

   int pad = 20;
   fprintf(fp, "\n");
   fprintf(fp, "===== HMM PROFILE ====================================\n");
   fprintf(fp, "%*s\t%s\n",  pad,  "NAME",       prof->name);
   fprintf(fp, "%*s\t%d\n",  pad,  "LENGTH",     prof->N);
   fprintf(fp, "%*s\t%d\n",  pad,  "ALLOC",      prof->Nalloc);
   fprintf(fp, "%*s\t%s\n",  pad,  "CONSENSUS",  prof->consensus);

   /* background model */
   fprintf(fp, "#%10s:\t", "FREQ");
   for (int j = 0; j < NUM_AMINO; j++)
      fprintf(fp, "%9.4f ", prof->bg_model->freq[j]);
   fprintf(fp, "\n");

   fprintf(fp, "#%10s:\t", "COMPO");
   for (int j = 0; j < NUM_AMINO; j++)
      fprintf(fp, "%9.4f ", prof->bg_model->compo[j]);
   fprintf(fp, "\n");

   fprintf(fp, "#%10s:\t", "INSERT");
   for (int j = 0; j < NUM_AMINO; j++)
      fprintf(fp, "%9.4f ", prof->bg_model->insert[j]);
   fprintf(fp, "\n");

   fprintf(fp, "#%10s:\t", "TRANS");
   for (int j = 0; j < NUM_TRANS_STATES; j++)
      fprintf(fp, "%9.4f ", prof->bg_model->trans[j]);
   fprintf(fp, "\n\n");

   /* position-specific probabilities */
   fprintf(fp, " %5s ", "");
   for (int j = 0; j < NUM_AMINO_PLUS_SPEC; j++) {
      fprintf(fp, "%9d ", j);
   }
      fprintf(fp, "\n");
   for (int i = 0; i < prof->N+1; i++)
   {
      /* line 1: match emissions */
      fprintf(fp, " %5d ", i);
      for (int j = 0; j < NUM_AMINO_PLUS_SPEC; j++)
         fprintf(fp, "%9.4f ", prof->hmm_model[i].match[j]);
      fprintf(fp, "\n");

      /* line 2: insert emissions */
      fprintf(fp, "%7s", "");
      for (int j = 0; j < NUM_AMINO_PLUS_SPEC; j++)
         fprintf(fp, "%9.4f ", prof->hmm_model[i].insert[j]);
      fprintf(fp, "\n");
      
      /* line 3: transition probabilities */
      fprintf(fp, "%7s", "");
      for (int j = 0; j < NUM_TRANS_STATES; j++)
         fprintf(fp, "%9.4f ", prof->hmm_model[i].trans[j]);
      fprintf(fp, "\n");
   }
   fprintf(fp, "\n");

   fprintf(fp, "#%15s:\n", "BACKGROUND");

   fprintf(fp, "#%15s:%7s", "LOG", "");
   for (int i = 0; i < NUM_AMINO; i++)
      fprintf(fp, "%9.4f ", BG_MODEL_log[i]);
   fprintf(fp, "\n");

   fprintf(fp, "#%15s:%7s", "ACTUAL", "");
   for (int i = 0; i < NUM_AMINO; i++)
      fprintf(fp, "%9.4f ", BG_MODEL[i]);
   fprintf(fp, "\n");

   fprintf(fp, "#%15s:\n", "SPECIAL");
   fprintf(fp, "%16s:\t%9.4f %9.4f\n", "E", prof->bg_model->spec[SP_E][SP_LOOP], prof->bg_model->spec[SP_E][SP_MOVE]);
   fprintf(fp, "%16s:\t%9.4f %9.4f\n", "N", prof->bg_model->spec[SP_N][SP_LOOP], prof->bg_model->spec[SP_N][SP_MOVE]);
   fprintf(fp, "%16s:\t%9.4f %9.4f\n", "C", prof->bg_model->spec[SP_C][SP_LOOP], prof->bg_model->spec[SP_C][SP_MOVE]);
   fprintf(fp, "%16s:\t%9.4f %9.4f\n", "J", prof->bg_model->spec[SP_J][SP_LOOP], prof->bg_model->spec[SP_J][SP_MOVE]);

   fprintf(fp, "//\n");
}

