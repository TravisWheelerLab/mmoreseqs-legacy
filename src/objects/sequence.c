/*******************************************************************************
 *  @file sequence.c
 *  @brief SEQUENCE Object
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
   seq->filename = NULL;
   seq->name     = NULL;
   seq->alph     = NULL;
   seq->seq      = NULL;

   return seq;
}

/* Destructor */
void SEQUENCE_Destroy(SEQUENCE *seq)
{
   free(seq->filename);
   free(seq->name);
   free(seq->alph);
   free(seq->seq);

   free(seq);
}

/* Set Sequence String to SEQUENCE and update length */
void SEQUENCE_Set_Seq(SEQUENCE* seq,
                      char*     seq_text)
{
   seq->N   = strlen(seq_text);
   seq->seq = malloc( sizeof(char) * (strlen(seq_text) + 1) );
   if (seq->seq == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc SEQ_TEXT for SEQUENCE.\n");
      exit(EXIT_FAILURE);
   }
   strcpy( seq->seq, seq_text );
}

/* Append Sequence String onto current SEQUENCE and update length */
void SEQUENCE_Append_Seq(SEQUENCE* seq,
                         char*     seq_text)
{
   seq->N  += strlen(seq_text);
   seq->seq = realloc( seq->seq, sizeof(char) * (seq->N + 1) );
   if (seq->seq == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc SEQ_TEXT for SEQUENCE.\n");
      exit(EXIT_FAILURE);
   }
   strcat( seq->seq, seq_text );
}

/* Set Textfield to SEQUENCE field */
void SEQUENCE_Set_Textfield(char** seq_field,
                            char*  text)
{
   *seq_field = realloc( *seq_field, sizeof(char) * (strlen(text) + 1) );
   if (*seq_field == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc TEXTFIELD for SEQUENCE.\n");
      exit(EXIT_FAILURE);
   }
   strcat( *seq_field, text );
}

/* Output SEQUENCE to FILE POINTER */
void SEQUENCE_Dump(SEQUENCE* seq,
                   FILE*     fp)
{
   /* check for bad pointer */
   if (fp == NULL) {
      fprintf(stderr, "ERROR: Bad FILE POINTER for printing SEQUENCE.\n");
      exit(EXIT_FAILURE);
   }

   /* space padding */
   const int pad = 10;

   fprintf(fp, "===== SEQUENCE =========================================\n");
   fprintf(fp, "\t%*s:\t%s\n", pad, "NAME",      seq->name);
   fprintf(fp, "\t%*s:\t%d\n", pad, "LENGTH",    seq->N);
   fprintf(fp, "\t%*s:\t%s\n", pad, "SEQUENCE",  seq->seq);
   fprintf(fp, "========================================================\n");
   fprintf(fp, "\n");
}


