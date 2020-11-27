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


/** FUNCTION:  SEQUENCE_Create()
 *  SYNOPSIS:  
 *
 *  RETURN:    Return pointer to new SEQUENCE object.
 */
SEQUENCE* 
SEQUENCE_Create() 
{
   SEQUENCE *seq = ERROR_malloc(sizeof(SEQUENCE));

   seq->full_N   = 0;
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


/** FUNCTION:  SEQUENCE_Destroy()
 *  SYNOPSIS:  
 *
 *  RETURN:    If successful, returns NULL pointer.
 */
SEQUENCE* 
SEQUENCE_Destroy( SEQUENCE *seq )
{
   if ( seq == NULL ) return seq;

   ERROR_free(seq->header);
   ERROR_free(seq->filename);
   ERROR_free(seq->name);
   ERROR_free(seq->alph);

   ERROR_free(seq->seq);
   ERROR_free(seq->dseq);

   ERROR_free(seq);
   
   return NULL;
}


/** FUNCTION:  SEQUENCE_Reuse()
 *  SYNOPSIS:  Reuse sequence by reinitializing all fields except <seq> field.
 *
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors. 
 */
void 
SEQUENCE_Reuse( SEQUENCE* seq )
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


/** FUNCTION:  SEQUENCE_Reuse()
 *  SYNOPSIS:  Set Sequence String <seq_text> to SEQUENCE <seq> and update length.
 *
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors. 
 */
void 
SEQUENCE_Set_Seq( SEQUENCE* seq,
                  char*     seq_text )
{
   seq->N = strlen(seq_text);
   /* resize string if necessary */
   if ( seq->N >= seq->Nalloc ) {
      SEQUENCE_Resize_Seq( seq, (seq->N + 1) );
   }  
   strcpy( seq->seq, seq_text );
}


/** FUNCTION:  SEQUENCE_Append_Seq()
 *  SYNOPSIS:  Append Sequence String <seq_text> onto current SEQUENCE <seq> and update length.
 *
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors. 
 */
void 
SEQUENCE_Append_Seq( SEQUENCE* seq,
                     char*     seq_text )
{
   seq->N  += strlen(seq_text);
   /* resize string if necessary (allocate twice the space currently needed) */
   if ( seq->N >= seq->Nalloc ) 
      SEQUENCE_Resize_Seq( seq, 2 * (seq->N + 1) );
   strcat( seq->seq, seq_text );
}


/** FUNCTION:  SEQUENCE_Resize_Seq()
 *  SYNOPSIS:  Reallocate space for SEQUENCE <seq>
 *
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors. 
 */
void 
SEQUENCE_Resize_Seq( SEQUENCE* seq,
                     int       size )
{
   seq->Nalloc       = size;
   seq->seq          = ERROR_realloc( seq->seq, sizeof(char) * seq->Nalloc );
   seq->full_seq     = seq->seq;
   /* make sure final char in string is terminal char */
   seq->seq[size-1]  = '\0';
}


/** FUNCTION:  SEQUENCE_Resize_Seq()
 *  SYNOPSIS:  Set Textfield <test> to SEQUENCE field <seq_field> (overwrites).
 *
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors. 
 */
void 
SEQUENCE_Set_Textfield( char** seq_field,
                        char*  text )
{
   *seq_field = ERROR_realloc( *seq_field, sizeof(char) * (strlen(text) + 1) );
   strcpy( *seq_field, text );
}


/** FUNCTION:  SEQUENCE_Digitize()
 *  SYNOPSIS:  Digitize text <seq> to create digital sequence <dseq>. 
 *
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors. 
 */
void 
SEQUENCE_Digitize( SEQUENCE* seq )
{
   seq->dseq = ERROR_realloc( seq->seq, sizeof(int) * seq->Nalloc );

   for ( int i = 0; i < seq->N; i++ ) {
      seq->dseq[i] = seq->full_seq[i] - 'A';
   }
}


/** FUNCTION:  SEQUENCE_SetSubseq()
 *  SYNOPSIS:  Set SEQUENCE <seq> to cover a subsequence <seq> to sequence <full_seq>.
 *             Subsequence covers range <q_beg, q_end>.
 *
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors. 
 */
void 
SEQUENCE_SetSubseq(  SEQUENCE*  seq, 
                     int        q_beg,
                     int        q_end )
{
   seq->seq = seq->full_seq + q_beg;
   seq->N   = strlen( seq->seq );
}


/** FUNCTION:  SEQUENCE_UnsetSubseq()
 *  SYNOPSIS:  Set SEQUENCE <seq> to cover a subsequence <seq> to sequence <full_seq>.
 *             Subsequence covers range <q_beg, q_end>.
 *
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors. 
 */
void 
SEQUENCE_UnsetSubseq( SEQUENCE*  seq )
{
   seq->seq = seq->full_seq;
   seq->N   = strlen( seq->seq );
}


/** FUNCTION:  SEQUENCE_Dump()
 *  SYNOPSIS:  Set SEQUENCE <seq> to file pointer <fp>.
 *
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors. 
 */
void 
SEQUENCE_Dump( SEQUENCE* seq,
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


