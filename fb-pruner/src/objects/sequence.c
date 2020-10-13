/*******************************************************************************
 *  FILE:      sequence.c
 *  PURPOSE:   SEQUENCE object
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
#include <time.h>

/* local imports */
#include "structs.h"
#include "../utilities/utilities.h"
#include "objects.h"

/* header */
#include "sequence.h"


/* Constructor */
SEQUENCE* SEQUENCE_Create() 
{
   SEQUENCE *seq = malloc(sizeof(SEQUENCE));
   if (seq == NULL) {
      perror("Error while malloc'ing SEQUENCE.\n");
      exit(EXIT_FAILURE);
   }

   seq->N        = 0;
   seq->Nalloc   = 0;

   seq->header   = NULL;
   seq->filename = NULL;
   seq->name     = NULL;
   seq->alph     = NULL;

   seq->seq      = NULL;
   seq->dseq     = NULL;

   SEQUENCE_Resize_Seq( seq, 256 );

   return seq;
}

/* Destructor */
void* SEQUENCE_Destroy( SEQUENCE *seq )
{
   if ( seq == NULL ) return seq;

   ERROR_free(seq->header);
   ERROR_free(seq->filename);
   ERROR_free(seq->name);
   ERROR_free(seq->alph);

   ERROR_free(seq->seq);
   ERROR_free(seq->dseq);

   ERROR_free(seq);
   seq = NULL;
   return seq;
}

/* Reuse sequence by reinitializing all fields except seq field */
void SEQUENCE_Reuse( SEQUENCE* seq )
{
   seq->N = 0;

   ERROR_free(seq->filename);
   ERROR_free(seq->name);
   ERROR_free(seq->alph);

   seq->filename = NULL;
   seq->name     = NULL;
   seq->alph     = NULL;

   /* making first char the terminal char sets string length to zero */
   seq->seq[0] = '\0';
}

/* Set Sequence String to SEQUENCE and update length */
void SEQUENCE_Set_Seq(  SEQUENCE* seq,
                        char*     seq_text )
{
   seq->N = strlen(seq_text);
   /* resize string if necessary */
   if ( seq->N >= seq->Nalloc ) 
      SEQUENCE_Resize_Seq( seq, (seq->N + 1) );
   strcpy( seq->seq, seq_text );
}

/* Append Sequence String onto current SEQUENCE and update length */
void SEQUENCE_Append_Seq(  SEQUENCE* seq,
                           char*     seq_text )
{
   seq->N  += strlen(seq_text);
   /* resize string if necessary (allocate twice the space currently needed) */
   if ( seq->N >= seq->Nalloc ) 
      SEQUENCE_Resize_Seq( seq, 2 * (seq->N + 1) );
   strcat( seq->seq, seq_text );
}

/* Reallocate space for SEQUENCE */
void SEQUENCE_Resize_Seq( SEQUENCE* seq,
                          int       size )
{
   seq->Nalloc = size;
   seq->seq    = realloc( seq->seq, sizeof(char) * seq->Nalloc );
   if (seq->seq == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc SEQ_TEXT for SEQUENCE.\n");
      exit(EXIT_FAILURE);
   }
   /* make sure final char in string is terminal char */
   seq->seq[size-1] = '\0';
}

/* Set Textfield to SEQUENCE field (overwrites) */
void SEQUENCE_Set_Textfield(  char** seq_field,
                              char*  text )
{
   *seq_field = ERROR_realloc( *seq_field, sizeof(char) * (strlen(text) + 1) );
   strcpy( *seq_field, text );
}

/* Create a digitized version of current sequence */
void SEQUENCE_Digitize( SEQUENCE* seq )
{
   
}

/* Output SEQUENCE to FILE POINTER */
void SEQUENCE_Dump(  SEQUENCE* seq,
                     FILE*     fp )
{
   /* check for bad pointer */
   if (fp == NULL) {
      fprintf(stderr, "ERROR: Bad FILE POINTER for printing SEQUENCE.\n");
      exit(EXIT_FAILURE);
   }

   /* space padding */
   const int pad = 10;

   fprintf(fp, "===== SEQUENCE =========================================\n");
   fprintf(fp, "\t%*s:\t%s\n", pad, "NAME",     seq->name);
   fprintf(fp, "\t%*s:\t%d\n", pad, "LENGTH",   seq->N);
   fprintf(fp, "\t%*s:\t%d\n", pad, "ALLOC",    seq->Nalloc);
   fprintf(fp, "\t%*s:\t%s\n", pad, "SEQUENCE", seq->seq);
   fprintf(fp, "========================================================\n");
   fprintf(fp, "\n");
}


