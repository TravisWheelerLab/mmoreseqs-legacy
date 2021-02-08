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
#include "../utilities/_utilities.h"
#include "_objects.h"

/* header */
#include "sequence.h"


/*! FUNCTION:  SEQUENCE_Create()
 *  SYNOPSIS:  Create <seq>, allocate memory and return pointer.
 */
SEQUENCE* 
SEQUENCE_Create() 
{
   SEQUENCE* seq  = NULL;
   seq            = ERROR_malloc(sizeof(SEQUENCE));
   int min_size   = 256;

   seq->full_N    = 0;
   seq->N         = 0;
   seq->Nalloc    = 0;

   seq->header    = NULL;
   seq->filename  = NULL;
   seq->name      = NULL;
   seq->alph      = NULL;

   seq->seq       = NULL;
   seq->dseq      = NULL;

   SEQUENCE_Resize_Seq( seq, min_size );

   return seq;
}

/*! FUNCTION:  SEQUENCE_Destroy()
 *  SYNOPSIS:  Destroy <seq>, free memory and return NULL pointer.
 */
SEQUENCE* 
SEQUENCE_Destroy( SEQUENCE *seq )
{
   if ( seq == NULL ) return seq;

   seq->header    = STR_Destroy( seq->header );
   seq->filename  = STR_Destroy( seq->filename );
   seq->name      = STR_Destroy( seq->name );
   seq->alph      = STR_Destroy( seq->alph );

   ERROR_free(seq->seq);
   ERROR_free(seq->dseq);

   ERROR_free(seq);
   
   return NULL;
}

/*! FUNCTION:  SEQUENCE_Reuse()
 *  SYNOPSIS:  Reuse <seq> by reinitializing all fields except <seq> field.
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors. 
 */
STATUS_FLAG
SEQUENCE_Reuse( SEQUENCE* seq )
{
   seq->N = 0;

   seq->header    = STR_Destroy( seq->header );
   seq->filename  = STR_Destroy( seq->filename );
   seq->name      = STR_Destroy( seq->name );
   seq->alph      = STR_Destroy( seq->alph );

   /* making first char the terminal char sets string length to zero */
   seq->seq[0] = '\0';

   return STATUS_SUCCESS;
}

/** FUNCTION:  SEQUENCE_GetSeq()
 *  SYNOPSIS:  Get Sequence String <seq>.
 */
STR 
SEQUENCE_GetSeq(  SEQUENCE*   seq )
{
   return seq->seq;
}

/** FUNCTION:  SEQUENCE_SetSeq()
 *  SYNOPSIS:  Set Sequence String <seq_text> to SEQUENCE <seq> and update length.
 */
void 
SEQUENCE_SetSeq(  SEQUENCE*   seq,
                  char*       seq_text )
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
SEQUENCE_Append_Seq( SEQUENCE*   seq,
                     char*       seq_text )
{
   if (seq_text == NULL) return;
   
   seq->N  += strlen(seq_text);
   /* resize string if necessary (allocate twice the space currently needed) */
   if ( seq->N >= seq->Nalloc ) {
      SEQUENCE_Resize_Seq( seq, 2 * (seq->N + 1) );
   }
   strcat( seq->seq, seq_text );
}


/** FUNCTION:  SEQUENCE_Resize_Seq()
 *  SYNOPSIS:  Reallocate space for SEQUENCE <seq>
 *
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors. 
 */
void 
SEQUENCE_Resize_Seq( SEQUENCE*   seq,
                     int         size )
{
   seq->Nalloc       = size;
   seq->seq          = ERROR_realloc( seq->seq, sizeof(char) * seq->Nalloc );
   seq->full_seq     = seq->seq;
   /* make sure final char in string is terminal char */
   seq->seq[size - 1]  = NULL_CHAR;
}


/** FUNCTION:  SEQUENCE_Resize_Seq()
 *  SYNOPSIS:  Set Textfield <test> to SEQUENCE field <seq_field> (overwrites).
 *
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors. 
 */
void 
SEQUENCE_SetTextfield( STR*  seq_field,
                        STR   text )
{
   *seq_field = STR_Set( *seq_field, text );
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


/** FUNCTION:  SEQUENCE_SetDomain()
 *  SYNOPSIS:  Set SEQUENCE <seq> to cover a subsequence <seq> to sequence <full_seq>.
 *             Subsequence covers domain range <q_beg, q_end>.
 *
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors. 
 */
void 
SEQUENCE_SetDomain( SEQUENCE*  seq, 
                     RANGE      Q_range )
{
   seq->seq = seq->full_seq + Q_range.beg;
   seq->N   = Q_range.end - Q_range.beg + 1;
}


/** FUNCTION:  SEQUENCE_UnsetSubseq()
 *  SYNOPSIS:  Set SEQUENCE <seq> to cover a subsequence <seq> to sequence <full_seq>.
 *             Subsequence covers range <q_beg, q_end>.
 *
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors. 
 */
void 
SEQUENCE_UnsetDomain( SEQUENCE*  seq )
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


