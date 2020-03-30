/*******************************************************************************
 *  FILE:      hmm_parser.h
 *  PURPOSE:   Parses .hmm files into HMM_PROFILE object
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

/* Parse .hmm file and builds a HMM_PROFILE object */
HMM_PROFILE* HMM_PROFILE_Parse( char*   _filename_,
                                long    offset )
{
   /* initialize HMM_PROFILE object */
   HMM_PROFILE*   prof           = HMM_PROFILE_Create( _filename_ );

   /* parser vars */
   FILE*          fp             = NULL;     
   int            line_count     = -1;       /* line counter of current line in file */
   char*          line_buf       = NULL;     /* pointer to start of buffered line */
   size_t         line_buf_size  = 0;        /* length of entire <line_buf> array */
   size_t         line_size      = 0;        /* length of current line in <line_buf> array */
   
   char*          line_ptr       = NULL;     /* pointer for splitting <line_buf> by delimiters */
   char*          token          = NULL;     /* pointer for iterating over <line_ptr> */
   char*          header         = NULL;     /* temp for current line header in <line_buf> */
   char*          field          = NULL;     /* temp for current line field in <line_buf> */
   float          value          = 0.0f;     /* temp for casting field to float */

   bool           inHeader       = true;     /* checks if are in the header or body of file */
   HMM_NODE*      curr_node      = NULL;     /* pointer to currently accessed node in the <hmm_profile> */

   /* open file */
   fp = fopen(_filename_, "r");
   /* check for file read error */
   if (fp == NULL)
   {
      char*  str = NULL;
      fprintf(stderr, "ERROR: Bad FILE POINTER for HMM PARSER => %s\n", _filename_ );
      exit(EXIT_FAILURE);
   }
   /* jump to start of desired hmm model */
   fseek(fp, offset, SEEK_SET);

   /* PARSE FILE HEADER FIELDS */
   /* read file line-by-line */
   while ( ( line_size = getline ( &line_buf, &line_buf_size, fp ) ), line_size != -1  && inHeader) /* read a line */
   {
      line_count++;

      /* remove newline from end of line */
      if (line_buf_size > 0 && line_buf[line_buf_size-1] == '\n') {
         line_buf[--line_buf_size] = '\0';
      }

      line_ptr = line_buf;
      header = strtok_r(line_ptr, " \n", &line_ptr);

      /* check which <header> data is being filled */
      if ( strcmp( header, "NAME" ) == 0 ) 
      {
         field = strtok_r(line_ptr, " \n", &line_ptr);
         HMM_PROFILE_Set_TextField( &prof->name, field );
      }
      else if ( strcmp( header, "ACC" ) == 0 ) 
      {
         field = strtok_r(line_ptr, " \n", &line_ptr);
         HMM_PROFILE_Set_TextField( &prof->acc, field );
      }
      else if ( strcmp( header, "DESC" ) == 0 ) 
      {
         field = strtok_r(line_ptr, " \n", &line_ptr);
         HMM_PROFILE_Set_TextField( &prof->acc, field );
      }
      else if ( strcmp( header, "LENG" ) == 0 ) 
      {
         field = strtok_r(line_ptr, " \n", &line_ptr);
         int num_nodes = atoi(field);
         HMM_PROFILE_Set_Model_Length(prof, num_nodes);
      }
      else if ( strcmp( header, "ALPH" ) == 0 ) 
      {
         field = strtok_r(line_ptr, " \n", &line_ptr);
         HMM_PROFILE_Set_Alphabet(prof, field);
      }
      else if ( strcmp( header, "RF" ) == 0 ) {}
      else if ( strcmp( header, "MM" ) == 0 ) {}

      else if ( strcmp( header, "STATS" ) == 0 ) 
      {
         field = strtok_r(line_ptr, " \n", &line_ptr); /* LOCAL */
         field = strtok_r(line_ptr, " \n", &line_ptr); /* distribution type */
         float param1 = atof( strtok_r(line_ptr, " \n", &line_ptr) );
         float param2 = atof( strtok_r(line_ptr, " \n", &line_ptr) );

         HMM_PROFILE_Set_Distribution_Params(prof, param1, param2, field);
      }
      /* COMPO is the background composition of the model */
      else if ( strcmp( header, "COMPO" ) == 0 ) 
      {
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

   /* PARSE NODE PROBABILITY VALUES */

   /* set first node to zero */
   curr_node = &prof->hmm_model[0];
   curr_node->match[0] = 1.0f;
   for ( int i = 1; i < NUM_AMINO; i++ ) {
      curr_node->match[0] = 0.0f;
   }

   /* read file node-by-node */
   for ( int j = 1; j < prof->N; j++ )
   {
      /* get first line in current node */
      line_size = getline ( &line_buf, &line_buf_size, fp );
      token = strtok ( line_buf, " \n" ); /* get first word (index) */
      // idx = atoi(token); /* index of line in sequence (unnecessary) */

      curr_node = &prof->hmm_model[j];  /* get current node */
      token = strtok(NULL, " \n");        /* get next word */

      /* LINE 1: Match Emission Line */
      for ( int j = 0; j < NUM_AMINO; j++)
      {
         /* check if valid float */
         if ( token[0] == '*' ) {
            curr_node->match[j] = 0.0f;
         } else {
            value = atof( token );
            value = negln2real( value );
            curr_node->match[j] = value;
         }

         token = strtok(NULL, " \n"); /* get next word */
      }

      /* LINE 2: Insert Emission Line */
      getline ( &line_buf, &line_buf_size, fp ); /* get next line */
      token = strtok ( line_buf, " \n" ); /* get first word */

      for ( int j = 0; j < NUM_AMINO; j++)
      {
         /* check if valid float */
         if ( token[0] == '*' ) {
            curr_node->insert[j] = 0.0f;
         } else {
            value = atof( token );
            value = negln2real( value );
            curr_node->insert[j] = value;
         }

         token = strtok(NULL, " \n"); /* get next word */
      }

      /* LINE 3: State Transition Line */
      getline ( &line_buf, &line_buf_size, fp ); /* get next line */
      token = strtok ( line_buf, " \n" ); /* get first word */

      for ( int j = 0; j < NUM_TRANS_STATES - 1; j++)
      {
         /* check if valid float */
         if ( token[0] == '*' ) {
            curr_node->trans[j] = 0.0f;
         } else {
            value = atof( token );
            value = negln2real( value );
            curr_node->trans[j] = value;
         }

         token = strtok(NULL, " \n"); /* get next word */
      }
   }

   /* initial node */
   curr_node->match[0] = 1.0f;
   for ( int i = 1; i < NUM_AMINO; i++ )
   {
      curr_node->match[i] = 0.0f;
   }

   /* free line buffer */
   free ( line_buf );
   /* close file */
   fclose ( fp );

   return prof;
}

/* configures HMM_PROFILE to account for background model */
/* modeled after HMMER p7_ProfileConfig() */
void HMM_PROFILE_Config(HMM_PROFILE* prof, 
                        int          mode)
{
   int      k     = 0; 
   int      x     = 0;
   float    Z     = 0.f;
   float*   occ   = NULL;

   occ = malloc( sizeof(float) * (prof->N + 1) );
   if (occ == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc for OCCUPANCY for PROFILE CONFIGURATION.\n");
      exit(EXIT_FAILURE);
   }

   HMM_PROFILE_CalcOccupancy(prof, occ);

   prof->isLocal = Test_IsLocal(mode);
   prof->isMultihit = Test_IsMulti(mode);

   /* set first node to zero */
   for (k = 0; k < NUM_AMINO; k++)
   {
      prof->hmm_model[0].match[k] = -INF;
      prof->hmm_model[0].insert[k] = -INF;
      prof->hmm_model[prof->N].insert[k] = -INF;
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
      for (k = 1; k <= prof->N; k++) {
         Z += occ[k] * (float) (prof->N - k + 1);
      }
      for (k = 1; k <= prof->N; k++) {
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

      for (k = 1; k < prof->N; k++)
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
      prof->hmm_model[prof->N].trans[x] = -INF;
   }

   /* Transition scores */
   for (k = 1; k < prof->N; k++) {
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

   for (k = 1; k < prof->N+1; k++) {
      /* each match score is log odds: log( match_prob / background_prob ) */
      for (x = 0; x < prof->alph_leng; x++)
      {
         prof->hmm_model[k].match[x] = log( (double) prof->hmm_model[k].match[x] / prof->bg_model->freq[x] );
      }
   }

   /* Insert Emission scores */
   for (x = 0; x < prof->alph_leng; x++) {
      prof->hmm_model[0].insert[x] = -INF; /* initial I-to-M impossible */

      for (k = 1; k < prof->N; k++) {
         prof->hmm_model[k].insert[x] = 0.0f;
      }

      prof->hmm_model[prof->N].insert[x] = -INF; /* initial I-to-M impossible */
   }

   /* Remaining specials, [NCJ][MOVE | LOOP], set by HMM_PROFILE_ReconfigLength() */
   HMM_PROFILE_ReconfigLength(prof, 100);
}

/* Calculates the Occupancy for the HMM_PROFILE */
void HMM_PROFILE_CalcOccupancy(HMM_PROFILE* prof, 
                               float*       occ)
{
   int      k;
   occ[0] = 0.;
   occ[1] = prof->hmm_model[0].trans[M2I] + prof->hmm_model[0].trans[M2M];
   for (k = 2; k <= prof->N; k++)
   {
      occ[k] = ( occ[k-1] * (prof->hmm_model[k-1].trans[M2I] + prof->hmm_model[k-1].trans[M2M]) ) + 
               ( (1 - occ[k-1]) * prof->hmm_model[k-1].trans[D2M] );
   }
}

/* Reconfigure the Length of the HMM_PROFILE */
void HMM_PROFILE_ReconfigLength(HMM_PROFILE*  prof, 
                                int           L)
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