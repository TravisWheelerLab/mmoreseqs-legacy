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
#include "structs.h"
#include "utilities.h"
#include "objects.h"

/* header */
#include "parsers.h"


/* Parse .hmm file and builds a HMM_PROFILE object */
void HMM_PROFILE_Parse( HMM_PROFILE*   prof,
                        char*          _filename_,
                        long           offset )
{
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

   const char*    delim          = " \t\n";  /* delimiters for parsing tokens (whitespace) */

   bool           inHeader       = true;     /* checks if are in the header or body of file */
   HMM_NODE*      curr_node      = NULL;     /* pointer to currently accessed node in the <hmm_profile> */

   const float    star_log       = NAN;      /* hmm profile stores values as logs, * represents zero b/c ln(0) is undefined */
   const float    star_real      = 0.f;      /* what to convert star(*) into in real space */

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

   /* reset the HMM_PROFILE fields for reuse */
   HMM_PROFILE_Reuse( prof );

   /* PARSE FILE HEADER FIELDS */
   /* read file line-by-line */
   while ( ( line_size = getline ( &line_buf, &line_buf_size, fp ) ), line_size != -1  && inHeader) /* read a line */
   {
      line_count++;

      line_ptr = line_buf;
      header = strtok(line_buf, delim);

      /* check which <header> data is being filled */
      if ( strcmp( header, "NAME" ) == 0 )
      {
         field = strtok(NULL, delim);
         // STRING_Replace(field, ' ', '_');
         HMM_PROFILE_Set_TextField( &prof->name, field );
      }
      else if ( strcmp( header, "ACC" ) == 0 )
      {
         field = strtok(NULL, delim);
         HMM_PROFILE_Set_TextField( &prof->acc, field );
      }
      else if ( strcmp( header, "DESC" ) == 0 )
      {
         field = strtok(NULL, delim);
         HMM_PROFILE_Set_TextField( &prof->acc, field );
      }
      else if ( strcmp( header, "LENG" ) == 0 )
      {
         field = strtok(NULL, delim);
         int num_nodes = atoi(field);
         HMM_PROFILE_Set_Model_Length(prof, num_nodes);
      }
      else if ( strcmp( header, "ALPH" ) == 0 )
      {
         field = strtok(NULL, delim);
         HMM_PROFILE_Set_Alphabet(prof, field);
      }
      else if ( strcmp( header, "RF" ) == 0 ) {}
      else if ( strcmp( header, "MM" ) == 0 ) {}

      else if ( strcmp( header, "STATS" ) == 0 )
      {
         field = strtok(NULL, delim); /* LOCAL */
         field = strtok(NULL, delim); /* distribution type */

         float param1 = atof( strtok(NULL, delim) );
         float param2 = atof( strtok(NULL, delim) );

         HMM_PROFILE_Set_Distribution_Params(prof, param1, param2, field);
      }
      /* COMPO is the background composition of the hmm model */
      else if ( strcmp( header, "COMPO" ) == 0 )
      {
         /* LINE 1: optional args (do nothing) */
         token = strtok(NULL, delim);
         for ( int j = 0; j < NUM_AMINO && token != NULL; j++)
         {
            value = atof( token );

            /* check if valid float */
            if ( token[0] != '*' ) {
               prof->bg_model->compo[j] = value;
            } else {
               prof->bg_model->compo[j] = star_log;
            }

            token = strtok(NULL, delim);
         }

         /* LINE 2: background emission probs */
         line_size = getline ( &line_buf, &line_buf_size, fp ); /* get next line */
         token = strtok ( line_buf, delim ); /* get first word */

         for ( int j = 0; j < NUM_AMINO && token != NULL; j++)
         {
            value = atof( token );

            /* check if valid float */
            if ( token[0] != '*' ) {
               prof->bg_model->insert[j] = value;
               prof->hmm_model[0].insert[j] = value;
            } else {
               prof->bg_model->insert[j] = star_log;
               prof->hmm_model[0].insert[j] = star_log;
            }
            /* get next word */
            token = strtok(NULL, delim);
         }

         /* LINE 3: trans probs (identical for all positions) */
         line_size = getline ( &line_buf, &line_buf_size, fp );  /* get next line */
         token = strtok ( line_buf, delim ); /* get first word */

         for ( int j = 0; j < NUM_TRANS_STATES && token != NULL; j++ )
         {
            value = atof( token );

            /* check if valid float */
            if ( token[0] != '*' ) {
               prof->bg_model->trans[j] = value;
               prof->hmm_model[0].trans[j] = value;
            } else {
               prof->bg_model->trans[j] = star_log;
               prof->hmm_model[0].trans[j] = star_log;
            }
            /* get next word */
            token = strtok(NULL, delim);
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
   for ( int j = 1; j < (prof->N + 1); j++ )
   {
      /* get first line in current node */
      line_size = getline ( &line_buf, &line_buf_size, fp );
      /* get first word (index) */
      token = strtok ( line_buf, delim );
      /* index of line in sequence (unnecessary) */
      // idx = atoi(token);

      /* get current node */
      curr_node = &prof->hmm_model[j];
      /* get next word */
      token = strtok(NULL, delim);

      /* LINE 1: Match Emission Line */
      for ( int j = 0; j < NUM_AMINO; j++)
      {
         /* check if valid float */
         if ( token[0] == '*' ) {
            curr_node->match[j] = star_log;
         } else {
            value = atof( token );
            curr_node->match[j] = value;
         }
         /* get next word */
         token = strtok(NULL, delim);
      }

      /* LINE 2: Insert Emission Line */
      line_size = getline ( &line_buf, &line_buf_size, fp ); /* get next line */
      token = strtok ( line_buf, delim ); /* get first word */

      for ( int j = 0; j < NUM_AMINO; j++)
      {
         /* check if valid float */
         if ( token[0] == '*' ) {
            curr_node->insert[j] = star_log;
         } else {
            value = atof( token );
            curr_node->insert[j] = value;
         }
         /* get next word */
         token = strtok(NULL, delim);
      }

      /* LINE 3: State Transition Line */
      line_size = getline ( &line_buf, &line_buf_size, fp ); /* get next line */
      token = strtok ( line_buf, delim ); /* get first word */

      for ( int j = 0; j < NUM_TRANS_STATES - 1; j++)
      {
         /* check if valid float */
         if ( token[0] == '*' ) {
            curr_node->trans[j] = star_log;
         } else {
            value = atof( token );
            curr_node->trans[j] = value;
         }
         /* get next word */
         token = strtok(NULL, delim);
      }
   }

   // /* final node */
   // curr_node = &prof->hmm_model[0];
   // curr_node->match[0] = 1.0f;
   // for ( int i = 1; i < NUM_AMINO; i++ )
   // {
   //    curr_node->match[i] = 0.0f;
   // }

   /* free line buffer */
   free ( line_buf );
   /* close file */
   fclose ( fp );

   prof->numberFormat = PROF_FORMAT_NEGLOG;
}

/* .hmm stores numbers in log space, but we need reals */
void HMM_PROFILE_Convert_NegLog_To_Real( HMM_PROFILE* prof )
{
   /* verify profile is in negative log space */
   if ( prof->numberFormat != PROF_FORMAT_NEGLOG ) return;

   HMM_NODE*   curr_node   = NULL;
   /* value representing star(*) in hmm file (log space) */
   const float star_log    = NAN;
   /* value representing star(*) in hmm file (real space) */
   const float star_real   = 0.f;
   float       value       = 0.f;

   /* background model */
   for ( int j = 0; j < NUM_AMINO; j++ )
   {
      /* composition */
      value = prof->bg_model->compo[j];
      if ( isnan(value) ) {
         prof->bg_model->compo[j] = star_real;
      }
      else {
         prof->bg_model->compo[j] = negln2real( value );
      }
      /* emission */
      value = prof->bg_model->insert[j];
      if ( isnan(value) ) {
         prof->bg_model->insert[j] = star_real;
      }
      else {
         prof->bg_model->insert[j] = negln2real( value );
      }
   }
   for ( int j = 0; j < NUM_TRANS_STATES - 1; j++ )
   {
      /* transition */
      value = prof->bg_model->trans[j];
      if ( isnan(value) ) {
         prof->bg_model->trans[j] = star_real;
      }
      else {
         prof->bg_model->trans[j] = negln2real( value );
      }
   }

   /* hmm model */
   for ( int i = 0; i < (prof->N + 1); i++ )
   {
      curr_node = &prof->hmm_model[i];
      /* match emission */
      for ( int j = 0; j < NUM_AMINO; j++ )
      {
         value = curr_node->match[j];
         if ( isnan(value) ) {
            curr_node->match[j] = star_real;
         }
         else {
            curr_node->match[j] = negln2real( value );
         }
      }
      /* insert emission */
      for ( int j = 0; j < NUM_AMINO; j++ )
      {
         value = curr_node->insert[j];
         if ( isnan(value) ) {
            curr_node->insert[j] = star_real;
         }
         else {
            curr_node->insert[j] = negln2real( value );
         }
      }
      /* transition state */
      for ( int j = 0; j < NUM_TRANS_STATES - 1; j++ )
      {
         value = curr_node->trans[j];
         if ( isnan(value) ) {
            curr_node->trans[j] = star_real;
         }
         else {
            curr_node->trans[j] = negln2real( value );
         }
      }
   }

   /* set first node to zero */
   curr_node = &prof->hmm_model[0];
   curr_node->match[0] = 1.0f;
   for ( int i = 1; i < NUM_AMINO; i++ ) {
      curr_node->match[0] = 0.0f;
   }

   /* final node */
   curr_node = &prof->hmm_model[0];
   curr_node->match[0] = 1.0f;
   for ( int i = 1; i < NUM_AMINO; i++ )
   {
      curr_node->match[i] = 0.0f;
   }

   prof->numberFormat = PROF_FORMAT_REAL;
}

/* configures HMM_PROFILE to account for background model */
/* function modeled after HMMER p7_ProfileConfig() */
void HMM_PROFILE_Config( HMM_PROFILE* prof,
                         int          mode )
{
   int      k     = 0;
   int      x     = 0;
   float    Z     = 0.f;
   float*   mocc  = NULL;

   mocc = malloc( sizeof(float) * (prof->N + 1) );
   if (mocc == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc for OCCUPANCY for PROFILE CONFIGURATION.\n");
      exit(EXIT_FAILURE);
   }

   /* if profile is not already in real space, do so now */
   if ( prof->numberFormat == PROF_FORMAT_NEGLOG )
      HMM_PROFILE_Convert_NegLog_To_Real( prof );

   prof->isLocal = Test_IsLocal(mode);
   prof->isMultihit = Test_IsMulti(mode);

   /* set first node to zero */
   for (k = 0; k < NUM_AMINO; k++)
   {
      prof->hmm_model[0].match[k] = -INF;
      prof->hmm_model[0].insert[k] = -INF;
      prof->hmm_model[prof->N].insert[k] = -INF;
   }
   // for (k = 0; k < (NUM_TRANS_STATES - 1); k++)
   // {
   //    prof->hmm_model[0].trans[k] = -INF;
   // }

   /* Entry Scores... */
   if (prof->isLocal)
   {
      /*
       * Local mode entry:  occ[k] /( \sum_i occ[i] * (M-i+1))
       * (Reduces to uniform 2/(M(M+1)) for occupancies of 1.0)
       */
      HMM_PROFILE_CalcOccupancy(prof, mocc, NULL);

      /* set B-to-M transition probs (local entry mode) */
      Z = 0.f;
      for (k = 1; k <= prof->N; k++)
         Z += mocc[k] * (float) (prof->N - k + 1);
      for (k = 1; k <= prof->N; k++)
         prof->hmm_model[k - 1].trans[B2M] = (float) log( (double) mocc[k] / (double) Z );

      free(mocc);
   }
   else /* glocal modes: left wing retraction(?) */
   {
      Z = log(prof->hmm_model[0].trans[M2D]);
      prof->hmm_model[0].trans[B2M] = log( 1.0 - prof->hmm_model[0].trans[M2D] );
      // printf("%.3f\t", prof->hmm_model[0].trans[B2M]);

      for (k = 1; k < prof->N; k++)
      {
         prof->hmm_model[k].trans[B2M] = log( 1.0 - prof->hmm_model[k].trans[D2M] );
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
      prof->num_J                         = 1.0f;
   }
   else
   {
      prof->bg_model->spec[SP_E][SP_MOVE] = 0.0f;
      prof->bg_model->spec[SP_E][SP_LOOP] = -INF;
      prof->num_J                         = 0.0f;
   }

   /* Initialize Transition scores */
   for (x = 0; x < NUM_TRANS_STATES; x++) {
      prof->hmm_model[prof->N].trans[x] = -INF;
   }

   /* Transition scores */
   for (k = 1; k < prof->N; k++) {
      for (x = 0; x < NUM_TRANS_STATES - 1; x++) {
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

   for (k = 1; k < prof->N + 1; k++) {
      /* each match score is log odds: log( match_prob / background_prob ) */
      for (x = 0; x < prof->alph_leng; x++)
      {
         prof->hmm_model[k].match[x] = (float) log( (double) prof->hmm_model[k].match[x] / (double) prof->bg_model->freq[x] );
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

   /* Compute special transition probabilities */
   /* Temporary configuration to be updated when length of query sequence is known */
   HMM_PROFILE_ReconfigLength(prof, 100);

   prof->numberFormat = PROF_FORMAT_LOGODDS;
}

/* Calculates the Occupancy for the HMM_PROFILE */
void HMM_PROFILE_CalcOccupancy(HMM_PROFILE*  prof,
                               float*        mocc,
                               float*        iocc )
{
   int k = 0;

   if (mocc != NULL)
   {
      mocc[0] = 0.;
      mocc[1] = prof->hmm_model[0].trans[M2I] + prof->hmm_model[0].trans[M2M];
      for (k = 2; k <= prof->N; k++) {
         mocc[k] = ( mocc[k - 1] * (prof->hmm_model[k - 1].trans[M2I] + prof->hmm_model[k - 1].trans[M2M]) ) +
                   ( (1.0 - mocc[k - 1]) * prof->hmm_model[k - 1].trans[D2M] );
      }
   }

   if (iocc != NULL)
   {
      iocc[0] = prof->hmm_model[0].trans[M2I] / prof->hmm_model[0].trans[I2M];
      for (k = 1; k <= prof->N; k++) {
         iocc[k] = mocc[k] * prof->hmm_model[k].trans[M2I] / prof->hmm_model[k].trans[I2M];
      }
   }
}

/* Configure the Length of the HMM_PROFILE based on the length of the sequence */
void HMM_PROFILE_ReconfigLength(HMM_PROFILE*  prof,
                                int           L)
{
   float ploop, pmove;
   float nj = (float) prof->num_J;

   /* Configure N,J,C transitions so they bear L/(2+nj) of the total unannotated sequence length L. */
   pmove = (2.0f + nj) / ((float) L + 2.0f + nj);  /* 2/(L+2) for sw; 3/(L+3) for fs */
   ploop = 1.0f - pmove;

   /* hardwire numbers from p7_ReconfigLength() */
   prof->bg_model->spec[SP_N][SP_LOOP] = log( ploop );
   prof->bg_model->spec[SP_C][SP_LOOP] = log( ploop );
   prof->bg_model->spec[SP_J][SP_LOOP] = log( ploop );

   prof->bg_model->spec[SP_N][SP_MOVE] = log( pmove );
   prof->bg_model->spec[SP_C][SP_MOVE] = log( pmove );
   prof->bg_model->spec[SP_J][SP_MOVE] = log( pmove );
}