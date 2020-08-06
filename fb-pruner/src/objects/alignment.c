/*******************************************************************************
 *  FILE:      alignment.c
 *  PURPOSE:   ALIGNMENT Object.
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
#include "alignment.h"

/* constructor */
ALIGNMENT* ALIGNMENT_Create()
{
   ALIGNMENT *aln       = NULL;
   const int min_size   = 256;

   aln = (ALIGNMENT*) malloc( sizeof(ALIGNMENT) );
   if (aln == NULL) {
      const char* obj_name = "ALIGNMENT";
      fprintf(stderr, "ERROR: Couldn't malloc %s: <%p>.\n", obj_name, aln);
      exit(EXIT_FAILURE);
   }

   aln->seq_beg   = VECTOR_INT_Create();
   aln->seq_end   = VECTOR_INT_Create();
   aln->sequence  = VECTOR_CHAR_Create();

   aln->tr_beg    = VECTOR_INT_Create();
   aln->tr_end    = VECTOR_INT_Create();
   aln->traces    = VECTOR_TRACE_Create();

   aln->Q      = -1;
   aln->T      = -1;

   ALIGNMENT_Resize(aln, min_size);

   return aln;
}

/* destructor */
void* ALIGNMENT_Destroy(ALIGNMENT* aln)
{
   if (aln == NULL) return aln;
   
   VECTOR_INT_Destroy( aln->seq_beg );
   VECTOR_INT_Destroy( aln->seq_beg );
   VECTOR_CHAR_Destroy( aln->sequence );

   VECTOR_INT_Destroy( aln->tr_beg );
   VECTOR_INT_Destroy( aln->tr_end );
   VECTOR_TRACE_Destroy( aln->traces );

   free( aln );
   aln = NULL;

   return aln;
}

/* reuse alignment by clearing traceback and setting dimensions */
void ALIGNMENT_Reuse(ALIGNMENT*  aln,
                     int         Q,
                     int         T )
{
   VECTOR_INT_Reuse( aln->seq_beg );
   VECTOR_INT_Reuse( aln->seq_end );
   VECTOR_CHAR_Reuse( aln->sequence );

   VECTOR_INT_Create( aln->tr_beg );
   VECTOR_INT_Create( aln->tr_end );
   VECTOR_TRACE_Create( aln->traces );
}

/* push trace onto end of alignment */
void ALIGNMENT_Pushback(ALIGNMENT* aln,
                        TRACE*     tr)
{
   /* if debugging, do edgechecks */
   #if DEBUG
      /* if normal state, check bounds of normal dp matrix */
      /* if special state, check bounds of special dp matrix */
   #endif

   VECTOR_TRACE_Pushback( aln->traces, *tr );
}

/* resize TRACE array in ALIGNMENT */
void ALIGNMENT_Resize( ALIGNMENT*   aln,
                       int          size )
{
   VECTOR_TRACE_Resize( aln->traces, size );
}

/* resize TRACE array in ALIGNMENT */
void ALIGNMENT_GrowTo( ALIGNMENT*   aln,
                       int          size )
{
   VECTOR_TRACE_GrowTo( aln->traces, size );
}

/* Empty ALIGNMENT Array */
void ALIGNMENT_Clear(ALIGNMENT* aln)
{
   VECTOR_TRACE_Resize( aln->traces, 0 );
}

/* create string of alignment */
void ALIGNMENT_Build_String(  ALIGNMENT*     aln,
                              SEQUENCE*      query,
                              HMM_PROFILE*   target )
{
   int      i;                         /* counter for traceback index */
   int      q_0, t_0;                  /* query and target index of traceback */
   int      f_mat, l_mat;              /* first and last match */
   char     q_char, t_char, a_char;    /* query, target, and alignment character */
   char*    q_seq;                     /* query sequence */
   char*    t_seq;                     /* target sequence */
   TRACE*   tr = NULL;                 /* traceback */
   bool     in_aln = false;            /* currently in an alignment */
   int      N;                         /* length of alignment */

   N = aln->traces->N;

   /* need consensus string if it is not already built */
   if ( target->consensus == NULL ) {
      HMM_PROFILE_Set_Consensus( target );
   }

   /* if alignment consensus has not been create, do so now */
   VECTOR_CHAR_Reuse( aln->sequence );
   VECTOR_INT_Reuse( aln->seq_beg );
   VECTOR_INT_Reuse( aln->seq_end );

   VECTOR_INT_Reuse( aln->tr_beg );
   VECTOR_INT_Reuse( aln->tr_end );

   /* pass through the alignment until we find our first begin state */
   tr = aln->traces->data;
   i = 0;

   while ( i < N )
   {
      /* find start of alignment */
      for (i; i < N; i++, tr++)
      {
         q_0 = tr->i;
         t_0 = tr->j;

         /* at beginning of alignment in traceback and sequence string */
         if ( tr->st == B_ST ) {
            VECTOR_INT_Pushback( aln->tr_beg, i );
            VECTOR_INT_Pushback( aln->seq_beg, aln->sequence->N );
            in_aln = true;
            break;
         }
      }

      /* get beginning points of alignment */
      tr++;
      q_seq = &(query->seq[tr->i]);
      t_seq = &(target->consensus[tr->j]);

      /* get consensus of each alignment */
      for (i; i < N; i++, tr++, q_seq++, t_seq++)
      {
         /* at end of alignment */
         if ( tr->st == E_ST ) {
            VECTOR_INT_Pushback( aln->tr_end, i );
            VECTOR_INT_Pushback( aln->seq_end, aln->sequence->N );
            VECTOR_CHAR_Pushback( aln->sequence, '\0' );
            in_aln = false;
            break;
         }

         q_char = *q_seq;
         t_char = *t_seq;

         /* if characters match, report character */
         if ( q_char == t_char ) {
            a_char = q_char;
            VECTOR_CHAR_Pushback( aln->sequence, a_char );
            continue;
         }

         /* if one upper and one lower case, report lower case character */
         q_char = tolower( q_char );
         t_char = tolower( t_char );
         if ( q_char == t_char ) {
            a_char = q_char;
            VECTOR_CHAR_Pushback( aln->sequence, a_char );
            continue;
         }

         /* if it is a match state, report positive score */
         if ( tr->st == M_ST ) {
            a_char = '+';
            VECTOR_CHAR_Pushback( aln->sequence, a_char );
            continue;
         }

         /* if it is a insert state, report '-' score */
         if ( tr->st == I_ST ) {
            a_char = '|';
            VECTOR_CHAR_Pushback( aln->sequence, a_char );
            continue;
         }

         /* if it is a delete state, report '|' score */
         if ( tr->st == D_ST ) {
            a_char = '-';
            VECTOR_CHAR_Pushback( aln->sequence, a_char );
            continue;
         }
      }

   }
}

/* outputs ALIGNMENT to FILE pointer */
void ALIGNMENT_Save( ALIGNMENT*  aln,
                     char*       _filename_ )
{
   FILE* fp = fopen(_filename_, "w");
   ALIGNMENT_Dump(aln, fp);
   fclose(fp);
}

/* outputs ALIGNMENT to FILE pointer */
void ALIGNMENT_Dump( ALIGNMENT*  aln,
                     FILE*       fp )
{
   int      N  = aln->traces->N;
   TRACE*   tr = aln->traces->data;

   /* test for bad file pointer */
   if (fp == NULL) {
      const char* obj_name = "ALIGNMENT";
      fprintf(stderr, "ERROR: Bad FILE pointer for printing %s: <%p>.\n", obj_name, aln);
      exit(EXIT_FAILURE);
      return;
   }
   
   fprintf(fp, "# ALIGNMENT (length=%d)\n", N );
   for (unsigned int i = 0; i < N; ++i)
   {
      fprintf(fp, "[%d](%s,%d,%d)\n", i, STATE_NAMES[tr[i].st], tr[i].i, tr[i].j);
   }
}
