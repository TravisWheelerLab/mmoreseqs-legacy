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

// /* Proteins in Alphabetical order */
// static char AA[] = { 
//    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
//    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
//    /* special characters */
//    '-', 
//    /* degen characters */
//    'X', 
//    /* special characters */
//    '*', '~'
// };

// /* Maps ASCII Letters to corresponding letters in AA[] */
// static int AA_REV[] = {
//    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
//    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
//    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
//    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
//    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
//    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
//    -1, -1, -1, -1, -1,  0, -1,  1,  2,  3, /* begin uppercase alphabet */
//     4,  5,  6,  7, -1,  8, 9,  10, 11, -1,
//    12, 13, 14, 15, 16, -1, 17, 18, 20, 19, /* X - unknown character */
//    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 
// };

/*! FUNCTION:  SEQUENCE_Create()
 *  SYNOPSIS:  Create <seq>, allocate memory and return pointer.
 */
SEQUENCE* 
SEQUENCE_Create() 
{
   SEQUENCE* seq        = NULL;
   seq                  = ERROR_malloc( sizeof(SEQUENCE) );
   int min_size         = 256;

   seq->full_N          = 0;
   seq->N               = 0;
   seq->Nalloc          = 0;

   seq->header          = NULL;
   seq->filename        = NULL;
   seq->name            = NULL;
   seq->alph            = NULL;

   seq->seq             = NULL;
   seq->dseq            = VECTOR_INT_Create();
   seq->vecseq          = VECTOR_CHAR_Create();
   seq->is_digitized    = false;

   SEQUENCE_Resize( seq, min_size );

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
   seq->vecseq    = VECTOR_CHAR_Destroy( seq->vecseq );
   seq->dseq      = VECTOR_INT_Destroy( seq->dseq );

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

   seq->header          = STR_Destroy( seq->header );
   seq->filename        = STR_Destroy( seq->filename );
   seq->name            = STR_Destroy( seq->name );
   seq->alph            = STR_Destroy( seq->alph );

   VECTOR_CHAR_Reuse( seq->vecseq );
   VECTOR_INT_Reuse( seq->dseq );
   seq->is_digitized    = false;

   /* making first char the terminal char sets string length to zero */
   seq->seq[0] = '\0';

   return STATUS_SUCCESS;
}

/*! FUNCTION:  SEQUENCE_GetSize()
 *  SYNOPSIS:  Gets the length of the sequence <seq>.
 */
size_t
SEQUENCE_GetSize( SEQUENCE* seq )
{
   size_t N;
   N = VECTOR_CHAR_GetSize( seq->vecseq );
   return N;
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
   /* old sequence type */
   seq->N = strlen(seq_text);
   /* resize string if necessary */
   if ( seq->N >= seq->Nalloc ) {
      SEQUENCE_Resize( seq, (seq->N + 1) );
   }  
   strcpy( seq->seq, seq_text );

   /* new sequence type */
   VECTOR_CHAR_Reuse( seq->vecseq );
   VECTOR_CHAR_Append( seq->vecseq, seq_text, strlen(seq_text) );
}

/** FUNCTION:  SEQUENCE_AppendSeq()
 *  SYNOPSIS:  Append Sequence String <seq_text> onto current SEQUENCE <seq> and update length.
 *
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors. 
 */
void 
SEQUENCE_AppendSeq(  SEQUENCE*   seq,
                     char*       seq_text )
{
   if (seq_text == NULL) return;
   
   /* old sequence type */
   seq->N  += strlen(seq_text);
   /* resize string if necessary (allocate twice the space currently needed) */
   if ( seq->N >= seq->Nalloc ) {
      SEQUENCE_Resize( seq, 2 * (seq->N + 1) );
   }
   strcat( seq->seq, seq_text );

   /* new sequence type */
   VECTOR_CHAR_Append( seq->vecseq, seq_text, strlen(seq_text) );
}


/** FUNCTION:  SEQUENCE_Resize()
 *  SYNOPSIS:  Reallocate space for SEQUENCE <seq>
 *
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors. 
 */
void 
SEQUENCE_Resize(  SEQUENCE*   seq,
                  int         size )
{
   /* old sequence type */
   seq->Nalloc       = size;
   seq->seq          = ERROR_realloc( seq->seq, sizeof(char) * seq->Nalloc );
   seq->full_seq     = seq->seq;
   /* make sure final char in string is terminal char */
   seq->seq[size - 1]  = NULL_CHAR;

   /* new sequence type */
   VECTOR_CHAR_GrowTo( seq->vecseq, (size_t)size );
}

/** FUNCTION:  SEQUENCE_GetCharAt()
 *  SYNOPSIS:  Get character <c> from sequence <seq> at <i>th position.
 */
char 
SEQUENCE_GetCharAt(  const SEQUENCE*   seq,
                     const int         i )
{
   char c;
   c = VECTOR_CHAR_Get( seq->vecseq, i );
   return c;
}

/** FUNCTION:  SEQUENCE_GetDigitAt()
 *  SYNOPSIS:  Get int <val> from digitized sequence at <i>th position.
 *             Caller must call _Digitize() before this.
 */
int 
SEQUENCE_GetDigitAt(    const SEQUENCE*   seq,
                        const int         i )
{
   char val;
   val = VECTOR_INT_Get( seq->dseq, i );
   return val;
}

/** FUNCTION:  SEQUENCE_Resize()
 *  SYNOPSIS:  Set Textfield <test> to SEQUENCE field <seq_field> (overwrites).
 *
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors. 
 */
void 
SEQUENCE_SetTextfield(  STR*  seq_field,
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
   if (seq->is_digitized == true) {
      return;
   }

   int N = SEQUENCE_GetSize( seq );
   VECTOR_INT_GrowTo( seq->dseq, N );
   for ( int i = 0; i < N; i++ ) {
      /* convert character to int value via lookup table */
      char  ch    = VECTOR_CHAR_Get( seq->vecseq, i );
      int   val   = AA_REV[ch];
      VECTOR_INT_Set( seq->dseq, i, val );
   }

   seq->is_digitized = true;
}


/** FUNCTION:  SEQUENCE_SetDomain()
 *  SYNOPSIS:  Set SEQUENCE <seq> to cover a subsequence <seq> to sequence <full_seq>.
 *             Subsequence covers domain range <q_beg, q_end>.
 *
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors. 
 */
void 
SEQUENCE_SetDomain(  SEQUENCE*  seq, 
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
      ERRORCHECK_exit(EXIT_FAILURE);
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


