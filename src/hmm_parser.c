/*******************************************************************************
 *  @file parser.c
 *  @brief Parses .hmm, .fasta, and .submat files
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/

// imports
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

// local imports
#include "structs.h"
#include "misc.h"
#include "hmm_parser.h"

/* parse .hmm file and build HMM_PROFILE object */
void hmmprofile_Create(HMM_PROFILE *prof, char *_filename_)
{
   printf("building profile object...\n");

   /* parser objects */
   char * token;
   char * header;
   char * field;
   float value;
   bool inHeader = true;

   /* line reader objects */
   char *line_buf = NULL;
   size_t line_buf_size = 0;
   int line_count = 0;
   ssize_t line_size;

   HMM_NODE *curr_node;

   /* open file */
   FILE *fp = fopen(_filename_, "r");
   if (fp == NULL)
   {
      char * str = NULL;
      perror( "Error while reading file." );
      exit(EXIT_FAILURE);
   }

   /* NOTE: May need to change to variable length string fields */

   /* parse header */
   while ( ( line_size = getline ( &line_buf, &line_buf_size, fp ) ), line_size != -1  && inHeader) /* read a line */
   {
      /* remove newline from end of line */
      if (line_buf_size > 0 && line_buf[line_buf_size-1] == '\n') {
         line_buf[--line_buf_size] = '\0';
      }

      /* delimit words on space and remove newline */
      token = strtok ( line_buf, " \n" );
      header = token;

      /* write the line */
      // printf ("[%d] %zd: %s \n", line_count, line_size, token );

      /* header - field */
      if ( strcmp( header, "NAME" ) == 0 ) {
         field = strtok( NULL, " \n" );
         prof->name = malloc( sizeof (char) * 256 );
         strcpy( prof->name, field );
      }
      else if ( strcmp( header, "ACC" ) == 0 ) {
         field = strtok( NULL, " \n" );
         prof->acc = malloc( sizeof (char) * 256 );
         strcpy( prof->acc, field );
      }
      else if ( strcmp( header, "DESC" ) == 0 ) {
         field = strtok( NULL, " \n");
         prof->desc = malloc( sizeof (char) * 256 );
         strcpy( prof->desc, field );
      }
      else if ( strcmp( header, "LENG" ) == 0 ) {
         field = strtok( NULL, " \n" );
         int value = atoi(field);
         prof->leng = value;

         /* allocate nodes */
         prof->hmm_model = malloc( sizeof (HMM_NODE) * (prof->leng + 1) );
      }
      else if ( strcmp( header, "ALPH" ) == 0 ) {
         field = strtok( NULL, " \n" );
         prof->alph = malloc( sizeof (char) * 256 );
         strcpy( prof->alph, field );

         int i = 0;
         while (field[i] != '\0') {
            field[i] = tolower(field[i]);
            i++;
         }

         if ( strcmp(field, "amino") == 0 ) {
            prof->alph_leng = NUM_AMINO;
         }
         else if ( strcmp(field, "dna") == 0 ) {
            prof->alph_leng = NUM_DNA;
         }
         else {
            char * str = NULL;
            perror( "Invalid data type: " );
            exit(EXIT_FAILURE);
         }
      }
      else if ( strcmp( header, "RF" ) == 0 ) {}
      else if ( strcmp( header, "MM" ) == 0 ) {}

      else if ( strcmp( header, "STATS" ) == 0 ) {
         field = strtok( NULL, " \n" ); /* LOCAL */
         field = strtok( NULL, " \n" ); /* distribution type */
         float param1 = atof( strtok( NULL, " \n" ) );
         float param2 = atof( strtok( NULL, " \n" ) );

         float *parPtr1 = NULL;
         float *parPtr2 = NULL;

         if ( strcmp( header, "MSV" ) == 0 ) {
            parPtr1 = &prof->msv_dist->param1;
            parPtr2 = &prof->msv_dist->param2;
         }
         else if ( strcmp( header, "VITERBI" ) == 0 ) {
            parPtr1 = &prof->viterbi_dist->param1;
            parPtr2 = &prof->viterbi_dist->param2;
         }
         else if ( strcmp( header, "FORWARD" ) == 0 ) {
            parPtr1 = &prof->forward_dist->param1;
            parPtr2 = &prof->forward_dist->param2;
         }

         parPtr1 = &param1;
         parPtr2 = &param2;
      }
      /* COMPO is the background composition of the model */
      else if ( strcmp( header, "COMPO" ) == 0 ) {

         /* malloc bg_model */
         prof->bg_model = malloc( sizeof (HMM_BG) );

         /* LINE 1: optional args (do nothing) */
         token = strtok ( NULL, " \n" );
         for ( int j = 0; j < NUM_AMINO && token != NULL; j++)
         {
            value = atof( token );
            value = negln2real( value );

            /* check if valid float */
            if ( token[0] != '*' ) {
               prof->bg_model->compo[j] = value;
            } else {
               prof->bg_model->compo[j] = 0.0f;
            }

            token = strtok ( NULL, " \n" );
         }

         /* LINE 2: background emission probs */
         getline ( &line_buf, &line_buf_size, fp ); /* get next line */
         token = strtok ( line_buf, " \n" ); /* get first word */

         for ( int j = 0; j < NUM_AMINO && token != NULL; j++)
         {
            value = atof( token );
            value = negln2real( value );

            /* check if valid float */
            if ( token[0] != '*' ) {
               prof->bg_model->insert[j] = value;
               prof->hmm_model[0].insert[j] = value;
            } else {
               prof->bg_model->insert[j] = 0.0f;
               prof->hmm_model[0].insert[j] = 0.0f;
            }

            token = strtok(NULL, " \n"); /* get next word */
         }

         /* LINE 3: trans probs (identical for all positions) */
         getline ( &line_buf, &line_buf_size, fp );  /* get next line */
         token = strtok ( line_buf, " \n" ); /* get first word */

         for ( int j = 0; j < NUM_TRANS_STATES && token != NULL; j++ )
         {
            value = atof( token );
            value = negln2real( value );

            /* check if valid float */
            if ( token[0] != '*' ) {
               prof->bg_model->trans[j] = value;
               prof->hmm_model[0].trans[j] = value;
            } else {
               prof->bg_model->trans[j] = 0.0f;
               prof->hmm_model[0].trans[j] = 0.0f;
            }

            token = strtok(NULL, " \n"); /* get next word */
         }

         inHeader = false;
      }

      line_count++;
      if (!inHeader) { break; }
   }

   /* malloc hmm_model (all nodes)  */
   int num_nodes = prof->leng;

   /* allocate current node */
   curr_node = &prof->hmm_model[0];
   curr_node->match[0] = 1.0f;

   /* parse hmm */
   int cnt = 1;
   int idx = 0;
   while ( ( line_size = getline ( &line_buf, &line_buf_size, fp ) ), line_size != -1 ) /* read a line */
   {
      token = strtok ( line_buf, " \n" ); /* get first word */
      idx = atoi(token); /* index of line in sequence */

      /* allocate current node */
      curr_node = &prof->hmm_model[cnt];

      token = strtok(NULL, " \n"); /* get next word */

      /* LINE 1: Match Emission Line */
      for ( int i = 0; i < NUM_AMINO && token != NULL; i++)
      {
         value = atof( token );
         value = negln2real( value );

         /* check if valid float */
         if ( token[0] == '*' ) {
            curr_node->match[i] = 0.0f;
         } else {
            curr_node->match[i] = value;
         }

         token = strtok(NULL, " \n"); /* get next word */
      }

      /* LINE 2: Insert Emission Line */
      getline ( &line_buf, &line_buf_size, fp ); /* get next line */
      token = strtok ( line_buf, " \n" ); /* get first word */

      for ( int i = 0; i < NUM_AMINO && token != NULL; i++)
      {
         value = atof( token );
         value = negln2real( value );

         /* check if valid float */
         if ( token[0] == '*' ) {
            curr_node->insert[i] = 0.0f;
         } else {
            curr_node->insert[i] = value;
         }

         token = strtok(NULL, " \n"); /* get next word */
      }

      /* LINE 3: State Transition Line */
      getline ( &line_buf, &line_buf_size, fp ); /* get next line */
      token = strtok ( line_buf, " \n" ); /* get first word */
      for ( int i = 0; i < NUM_TRANS_STATES && token != NULL; i++)
      {
         value = atof( token );
         value = negln2real( value );

         /* check if valid float */
         if ( token[0] == '*' ) {
            curr_node->trans[i] = 0.0f;
         } else {
            curr_node->trans[i] = value;
         }

         token = strtok(NULL, " \n"); /* get next word */
      }
      cnt++;
   }

   /* initial node */
   curr_node->match[0] = 1.0f;
   for ( int i = 1; i < NUM_AMINO; i++)
   {
      curr_node->match[i] = 0.0f;
   }

   /* free line buffer */
   free ( line_buf );
   line_buf = NULL;

   /* close file */
   fclose ( fp );
}

/* display HMM_PROFILE to console */
void hmmprofile_Save(HMM_PROFILE *prof, char *_filename_)
{
   int i,j;
   float best_val, new_val;
   char best_amino;
   HMM_NODE curr_node;

   FILE *fp = fopen(_filename_, "w");

   fprintf(fp, "\n");
   fprintf(fp, "===== HMM PROFILE ====================================\n");
   fprintf(fp, "\t         NAME:\t%s\n", prof->name);
   fprintf(fp, "\t       LENGTH:\t%d\n", prof->leng);
   char *seq = (char *)malloc( sizeof(char) * prof->leng );

   /* print consensus sequence */
   fprintf(fp, "\tCONSENSUS SEQ:\t");
   for (int i = 1; i < prof->leng+1; i++)
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
   // fprintf(fp, "\tHMM:\n\t");
   // for (int i = 0; i < NUM_AMINO; i++ {
   //    fprintf(fp, "%d\t", i);
   // }
   // fprintf(fp, "\n");
   fprintf(fp, "=== TRANSITION PROBS ===\n");
   for (int i = 0; i < prof->leng+1; i++)
   {
      fprintf(fp, "%d", i);
      for (int j = 0; j < NUM_TRANS_STATES; j++)
      {
         fprintf(fp, "\t%.4f", prof->hmm_model[i].trans[j]);
      }
      fprintf(fp, "\n");
   }
   fprintf(fp, "=== MATCH EMISSION PROBS ===\n");
   for (int i = 0; i < prof->leng+1; i++)
   {
      fprintf(fp, "%d", i);
      for (int j = 0; j < NUM_AMINO; j++)
      {
         fprintf(fp, "\t%.4f", prof->hmm_model[i].match[j]);
      }
      fprintf(fp, "\n");
   }
   fprintf(fp, "=== INSERT EMISSION PROBS ===\n");
   for (int i = 0; i < prof->leng+1; i++)
   {
      fprintf(fp, "%d", i);
      for (int j = 0; j < NUM_AMINO; j++)
      {
         fprintf(fp, "\t%.4f", prof->hmm_model[i].insert[j]);
      }
      fprintf(fp, "\n");
   }

   /* */
   // for (int i = 0; i < prof->leng+1; i++)
   // {
   //    fprintf(fp, "\t%d", i);
   //    for (int j = 0; j < NUM_AMINO; j++)
   //    {
   //       fprintf(fp, "\t%.4f", prof->hmm_model[i].match[j]);
   //    }
   //    fprintf(fp, "\n\t");
   //    for (int j = 0; j < NUM_AMINO; j++)
   //    {
   //       fprintf(fp, "\t%.4f", prof->hmm_model[i].insert[j]);
   //    }
   //    fprintf(fp, "\n\t");
   //    for (int j = 0; j < NUM_TRANS_STATES; j++)
   //    {
   //       fprintf(fp, "\t%.4f", prof->hmm_model[i].trans[j]);
   //    }
   //    fprintf(fp, "\n");
   // }
   // fprintf(fp, "\n");

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

   fclose(fp);
}

/* display HMM_PROFILE to console */
void hmmprofile_Display(HMM_PROFILE *prof)
{
   int i,j;
   float best_val, new_val;
   char best_amino;
   HMM_NODE curr_node;

   printf("\n");
   printf("===== HMM PROFILE ====================================\n");
   printf("\t         NAME:\t%s\n", prof->name);
   printf("\t       LENGTH:\t%d\n", prof->leng);
   char *seq = (char *)malloc( sizeof(char) * prof->leng );

   /* print consensus sequence */
   printf("\tCONSENSUS SEQ:\t");
   for (int i = 1; i < prof->leng+1; i++)
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
      printf("%c", best_amino);
   }
   printf("\n");

   printf("\t        COMPO:\n");
   printf("FREQ:\t");
   for (int j = 0; j < NUM_AMINO; j++)
   {
      printf("\t%.3f", prof->bg_model->freq[j]);
   }
   printf("\nCOMPO:\t");
   for (int j = 0; j < NUM_AMINO; j++)
   {
      printf("\t%.3f", prof->bg_model->compo[j]);
   }
   printf("\nINSERT:\t");
   for (int j = 0; j < NUM_AMINO; j++)
   {
      printf("\t%.3f", prof->bg_model->insert[j]);
   }
   printf("\nTRANS:\t");
   for (int j = 0; j < NUM_TRANS_STATES; j++)
   {
      printf("\t%.3f", prof->bg_model->trans[j]);
   }
   printf("\n\n");

   printf("\tHMM:\n\t\t");
   for (int i = 0; i < NUM_AMINO; i++)
   {
      printf("%d\t", i);
   }
   printf("\n");
   for (int i = 0; i < prof->leng+1; i++)
   {
      printf("\t%d", i);
      for (int j = 0; j < NUM_AMINO; j++)
      {
         printf("\t%.3f", prof->hmm_model[i].match[j]);
      }
      printf("\n\t");
      for (int j = 0; j < NUM_AMINO; j++)
      {
         printf("\t%.3f", prof->hmm_model[i].insert[j]);
      }
      printf("\n\t");
      for (int j = 0; j < NUM_TRANS_STATES; j++)
      {
         printf("\t%.3f", prof->hmm_model[i].trans[j]);
      }
      printf("\n\n");
   }
   printf("\n\n");

   printf("BACKGROUND:\n");
   printf("LOG:\t\t");
   for (int i = 0; i < NUM_AMINO; i++)
   {
      printf("%.3f\t", BG_MODEL_log[i]);
   }
   printf("\n");
   printf("ACTUAL:\t\t");
   for (int i = 0; i < NUM_AMINO; i++)
   {
      printf("%.3f\t", BG_MODEL[i]);
   }
   printf("\n\n");

   printf("SPECIAL:\n");
   printf("\tE:\t%.3f\t%.3f\n", prof->bg_model->spec[SP_E][SP_LOOP], prof->bg_model->spec[SP_E][SP_MOVE]);
   printf("\tN:\t%.3f\t%.3f\n", prof->bg_model->spec[SP_N][SP_LOOP], prof->bg_model->spec[SP_N][SP_MOVE]);
   printf("\tC:\t%.3f\t%.3f\n", prof->bg_model->spec[SP_C][SP_LOOP], prof->bg_model->spec[SP_C][SP_MOVE]);
   printf("\tJ:\t%.3f\t%.3f\n", prof->bg_model->spec[SP_J][SP_LOOP], prof->bg_model->spec[SP_J][SP_MOVE]);

   printf("\n======================================================\n\n");
}

/* configures HMM_PROFILE to account for background model */
/* modeled after HMMER p7_ProfileConfig() */
void hmmprofile_Config(HMM_PROFILE *prof, int mode)
{
   int k, x;
   float Z = 0.0;
   float *occ = malloc( sizeof(float) * (prof->leng + 1) );

   hmmprofile_CalcOccupancy(prof, occ);

   prof->isLocal = Test_IsLocal(mode);
   prof->isMultihit = Test_IsMulti(mode);

   /* set first node to zero */
   for (k = 0; k < NUM_AMINO; k++)
   {
      prof->hmm_model[0].match[k] = -INF;
      prof->hmm_model[0].insert[k] = -INF;
      prof->hmm_model[prof->leng].insert[k] = -INF;
   }
   for (k = 0; k < NUM_TRANS_STATES; k++)
   {
      prof->hmm_model[0].trans[k] = -INF;
   }

   /* Entry Scores... */
   // printf("B2M:\t");
   if (prof->isLocal)
   {
      /* set B-to-M transition probs (local entry mode) */
      for (k = 1; k <= prof->leng; k++) {
         Z += occ[k] * (float) (prof->leng - k + 1);
      }
      for (k = 1; k <= prof->leng; k++) {
         prof->hmm_model[k-1].trans[B2M] = log( occ[k] / Z );
         // printf("%.3f\t", prof->hmm_model[k-1].trans[B2M]);
      }
      free(occ);
   }
   else /* glocal modes: left wing retraction(?) */
   {
      Z = log(prof->hmm_model[0].trans[M2D]);
      prof->hmm_model[0].trans[B2M] = log(1.0 - prof->hmm_model[0].trans[M2D]);
      // printf("%.3f\t", prof->hmm_model[0].trans[B2M]);

      for (k = 1; k < prof->leng; k++)
      {
         prof->hmm_model[k].trans[B2M] = log(1.0 - prof->hmm_model[k].trans[D2M]);
         Z += log(prof->hmm_model[k].trans[D2D]);
         // printf("%.3f\t", prof->hmm_model[0].trans[B2M]);
      }
   }

   /* E state loop/move probabilities: nonzero for MOVE allows loops/multihits
    * N,C,J transitions are set later by length config
    */
   if (prof->isMultihit)
   {
      prof->bg_model->spec[SP_E][SP_MOVE] = -CONST_LOG2;
      prof->bg_model->spec[SP_E][SP_LOOP] = -CONST_LOG2;
      prof->bg_model->num_J = 1.0f;
   }
   else 
   {
      prof->bg_model->spec[SP_E][SP_MOVE] = 0.0f;
      prof->bg_model->spec[SP_E][SP_LOOP] = -INF;
      prof->bg_model->num_J = 0.0f;
   }

   /* Initialize Transition scores */
   for (x = 0; x < NUM_TRANS_STATES; x++) {
      prof->hmm_model[prof->leng].trans[x] = -INF;
   }

   /* Transition scores */
   for (k = 1; k < prof->leng; k++) {
      for (x = 0; x < NUM_TRANS_STATES-1; x++) {
         prof->hmm_model[k].trans[x] = log( prof->hmm_model[k].trans[x] );
      }
   }

   /* Set Background scores */
   for (x = 0; x < prof->alph_leng; x++) {
      prof->bg_model->freq[x] = BG_MODEL[x];
   }

   /* Match Emission scores */
   for (x = 0; x < prof->alph_leng; x++) {
      // prof->hmm_model[0].match[x] = -INF;
   }

   for (k = 1; k < prof->leng+1; k++) {
      /* each match score is log odds: log( match_prob / background_prob ) */
      for (x = 0; x < prof->alph_leng; x++)
      {
         prof->hmm_model[k].match[x] = log( (double) prof->hmm_model[k].match[x] / prof->bg_model->freq[x] );
      }
   }

   /* Insert Emission scores */
   for (x = 0; x < prof->alph_leng; x++) {
      prof->hmm_model[0].insert[x] = -INF; /* initial I-to-M impossible */

      for (k = 1; k < prof->leng; k++) {
         prof->hmm_model[k].insert[x] = 0.0f;
      }

      prof->hmm_model[prof->leng].insert[x] = -INF; /* initial I-to-M impossible */
   }

   /* Remaining specials, [NCJ][MOVE | LOOP], set by hmmprofile_ReconfigLength() */
   hmmprofile_ReconfigLength(prof, 100);
}

/* */
void hmmprofile_CalcOccupancy(HMM_PROFILE *prof, float *occ)
{
   int k;

   occ[0] = 0.;
   occ[1] = prof->hmm_model[0].trans[M2I] + prof->hmm_model[0].trans[M2M];
   for (k = 2; k <= prof->leng; k++)
   {
      occ[k] = ( occ[k-1] * (prof->hmm_model[k-1].trans[M2I] + prof->hmm_model[k-1].trans[M2M]) ) + 
               ( (1 - occ[k-1]) * prof->hmm_model[k-1].trans[D2M] );
   }
}

/* */
void hmmprofile_ReconfigLength(HMM_PROFILE *prof, int L)
{
   float ploop, pmove;
   float nj = (float) prof->bg_model->num_J;

   /* Configure N,J,C transitions so they bear L/(2+nj) of the total unannotated sequence length L. */
   pmove = (2.0f + nj) / ((float) L + 2.0f + nj);  /* 2/(L+2) for sw; 3/(L+3) for fs */
   ploop = 1.0f - pmove;

   /* hardwire numbers from p7_ReconfigLength() */
   // ploop = -0.02956; 
   // pmove = -3.53612;

   prof->bg_model->spec[SP_N][SP_LOOP] = log( ploop );
   prof->bg_model->spec[SP_C][SP_LOOP] = log( ploop );
   prof->bg_model->spec[SP_J][SP_LOOP] = log( ploop );

   prof->bg_model->spec[SP_N][SP_MOVE] = log( pmove );
   prof->bg_model->spec[SP_C][SP_MOVE] = log( pmove );
   prof->bg_model->spec[SP_J][SP_MOVE] = log( pmove );
}

/* parse .fasta file and build SEQ object */
void seq_Create(SEQ *seq, char *_filename_)
{
   /* line reader objects */
   char *line_buf = NULL;
   size_t line_buf_size = 0;
   int line_count = 0;
   ssize_t line_size;

   int num_seqs = 0;
   int seq_len;

   /* open file */
   FILE *fp;
   fp = fopen(_filename_, "r");

   if (fp == NULL)
   {
      char * str = NULL;
      perror( sprintf( str, "Error while reading file: ") );
      exit(EXIT_FAILURE);
   }

   while ( ( line_size = getline ( &line_buf, &line_buf_size, fp ) ), line_size != -1 ) /* read a line */
   {
      /* remove newline from end of line */
      if( line_buf[line_size-1] == '\n' ) {
         line_buf[line_size-1] = '\0';
         line_size -= 1;
      }

      // capitalize all dna
      for (int i=0; i<line_size; i++) {
         line_buf[i] = toupper(line_buf[i]);
      }

      /* write line to console */
      // printf ("[%d] %zd: %s \n", line_count, line_size, line_buf );

      /* check if line is a header */
      if (line_buf[0] == '>')
      {
         num_seqs++;

         /* if onto next seq, reset counters */
         if (num_seqs > 1) {
            seq->leng = seq_len;
            return;
         }

         /* add header name  */
         seq->name = (char *)malloc( sizeof(char) * line_size );
         strcpy( seq->name, line_buf+1 );

         /* reset counters */
         seq_len = 0;
         seq->seq = (char *)malloc( sizeof(char) );
      }
      else /* otherwise, append to current sequence */
      {
         /* resize sequence to new length and concat new line onto end of sequence */

         seq_len += line_size;
         seq->seq = (char *)realloc( seq->seq, sizeof(char) * seq_len);
         strcat( seq->seq, line_buf );
      }

      line_count++;
   }

   /* last seq data */
   seq->leng = seq_len;

   /* free line buffer */
   free ( line_buf );
   line_buf = NULL;

   /* close file */
   fclose ( fp );
}

/* display FASTA to console */
void seq_Display(SEQ *seq)
{
   printf("\n");
   printf("===== SEQ ============================================\n");
   printf("\t    NAME:\t%s\n", seq->name);
   printf("\t  LENGTH:\t%d\n", seq->leng);
   printf("\tSEQUENCE:\t%s\n", seq->seq);
   printf("======================================================\n\n");
}

/* parse .submat file and build SUBMAT object */
void submat_Create(SUBMAT *submat, char *_filename_)
{
   /* allocate memory for substitution matrix scores */
   submat->scores = (float *)malloc( sizeof(float) * SUBMAT_SIZE );

   /* line reader objects */
   char *line_buf = NULL;
   size_t line_buf_size = 0;
   int line_count = 0;
   ssize_t line_size;

   /* line parser objects */
	char delim[] = "\t";
	char *parser;
	int x,y,key;
   float score;
	char a,b;
	x = 0; y = 0;

   /* open file */
   FILE *fp;
   fp = fopen(_filename_, "r");

	if (fp != NULL)
	{
      /* read a line */
		while ( ( line_size = getline ( &line_buf, &line_buf_size, fp ) ), line_size != -1 )
		{
			b = AA[y];
			x = 0; /* reset x at start of each line */

         /* write the line */
         printf ("[%d] %zd: %s \n", line_count, line_size, line_buf );

			if (line_buf[0] != '#') /* if not a comment, parse line */
			{
				bool isNull = (parser != NULL);
				parser = strtok(line_buf, delim);

				// next row value is the first element on line
				if (parser != NULL) {
					b = parser[0];
					parser = strtok(NULL, delim);
				}
				while (parser != NULL) {
					// get current table character
					a = AA[x];
					// cast digit to integer (if not a digit, output is '0')
					score = atof(parser);

					// check if cast score is valid
					if (!(score == 0 && parser[0] != '0')) {
						// map protein pair to int (both directions)
						key = submat_Keymap(a,b);
						submat->scores[key] = score;
						key = submat_Keymap(b,a);
						submat->scores[key] = score;
						x += 1;
					}
					// split on next tab
					parser = strtok(NULL, delim);
				}
			}
			y += 1;
		}
	}
	else
	{
		perror ( _filename_ );
	}
	fclose(fp);
}

/* find data location */
int submat_Keymap(char x, char y)
{
   int X = x - 'A';
	int Y = y - 'A';

	return (X * ALPHA_MAX) + Y;
}

/* get score from SUBMAT */
float submat_Get(SUBMAT *submat, char q, char t)
{
   int key = submat_Keymap(q,t);
	return submat->scores[key];
}

/* display SUBMAT to console */
void submat_Display(SUBMAT *submat)
{
   printf("\n");
   printf("===== SUBMAT =========================================\n");

   char a,b;
   float score;
   int i,j = 0;

   printf("\t");
   for (i = 0; i < NUM_AMINO; i++)
   {
      printf("[%c]\t", AA[i]);
   }
   printf("\n");

   for (i = 0; i < NUM_AMINO; i++)
   {
      a = AA[i];
      printf("[%c]\t", a);
      for (j = 0; j < NUM_AMINO; j++)
      {
         b = AA[j];

         score = submat_Get(submat, a, b);
         printf("%.3f\t", score);
      }
      printf("\n");
   }

   printf("======================================================\n\n");
}

/* print RESULTS object to console */
void results_Display(RESULTS *results)
{
   HIT *hit;

   printf("\n");
   printf("===== RESULTS ========================================\n");
   printf("NUM_HITS: %d\n", results->N);
   for (int i = 0; i < results->N; i++)
   {
      hit = &results->hits[i];
      printf("> HIT[%d]:\n", i);
      printf("\t    QUERY:\t%s\n", hit->query->name);
      printf("\t   TARGET:\t%s\n", hit->target->name);
      printf("\tLOG-SCORE:\t%f\n", hit->score);
      printf("\t    SCORE:\t%f\n", negln2real(hit->score));
   }
   printf("======================================================\n\n");
}

/* convert negative natural logs to real probabilities */
float negln2real(float negln_prob)
{
   float real_prob = 1.0f / ( exp( negln_prob ) );
   return real_prob;
}

/* convert real probabilities to negative natural logs */
float real2negln(float real_prob)
{
   float negln_prob = -1.0f * log(real_prob);
   return negln_prob;
}
