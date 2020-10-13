/*******************************************************************************
 *  FILE:      alignment.c
 *  PURPOSE:   ALIGNMENT Object.
 *             Used for storing Viterbi or Forward/Backward Alignments.
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
#include "../utilities/utilities.h"
#include "objects.h"

/* header */
#include "alignment.h"

/*
 *  FUNCTION:  ALIGNMENT_Create()
 *  SYNOPSIS:
 */
ALIGNMENT* ALIGNMENT_Create()
{
   ALIGNMENT *aln       = NULL;
   const int min_size   = 256;

   aln = (ALIGNMENT*) ERROR_malloc( sizeof(ALIGNMENT) );

   aln->seq_beg      = VECTOR_INT_Create();
   aln->seq_end      = VECTOR_INT_Create();
   aln->sequence     = VECTOR_CHAR_Create();

   aln->tr_beg       = VECTOR_INT_Create();
   aln->tr_end       = VECTOR_INT_Create();
   aln->traces       = VECTOR_TRACE_Create();

   aln->Q            = -1;
   aln->T            = -1;
   aln->aln_len      = -1;

   aln->perc_id      = 0.0;
   aln->num_matches  = 0;
   aln->num_gaps     = 0;
   aln->num_misses   = 0;

   aln->cigar_aln    = NULL;
   aln->target_aln   = NULL;
   aln->center_aln   = NULL;
   aln->query_aln    = NULL;
   aln->state_aln    = NULL;

   ALIGNMENT_Resize(aln, min_size);

   return aln;
}

/*
 *  FUNCTION:  RESULTS_Destroy()
 *  SYNOPSIS:
 */
void* ALIGNMENT_Destroy(   ALIGNMENT* aln )
{
   if (aln == NULL) return aln;
   
   VECTOR_INT_Destroy( aln->seq_beg );
   VECTOR_INT_Destroy( aln->seq_end );
   VECTOR_CHAR_Destroy( aln->sequence );

   VECTOR_INT_Destroy( aln->tr_beg );
   VECTOR_INT_Destroy( aln->tr_end );
   VECTOR_TRACE_Destroy( aln->traces );

   ERROR_free( aln->cigar_aln );
   ERROR_free( aln->target_aln );
   ERROR_free( aln->center_aln );
   ERROR_free( aln->query_aln );
   ERROR_free( aln->state_aln );

   ERROR_free( aln );
   aln = NULL;

   return aln;
}

/*
 *  FUNCTION:  RESULTS_Reuse()
 *  SYNOPSIS:  Reuse alignment by clearing traceback and setting dimensions.
 */
void ALIGNMENT_Reuse(   ALIGNMENT*  aln,
                        int         Q,
                        int         T )
{
   VECTOR_INT_Reuse( aln->seq_beg );
   VECTOR_INT_Reuse( aln->seq_end );
   VECTOR_CHAR_Reuse( aln->sequence );

   VECTOR_INT_Reuse( aln->tr_beg );
   VECTOR_INT_Reuse( aln->tr_end );
   VECTOR_TRACE_Reuse( aln->traces );

   /* TODO: reuse string? */
   aln->cigar_aln = ERROR_free( aln->cigar_aln );
   aln->target_aln = ERROR_free( aln->target_aln );
   aln->center_aln = ERROR_free( aln->center_aln );
   aln->query_aln = ERROR_free( aln->query_aln );
   aln->state_aln = ERROR_free( aln->state_aln );
}

/* push trace onto end of alignment */
/*
 *  FUNCTION:  RESULTS_Reuse()
 *  SYNOPSIS:  Reuse alignment by clearing traceback and setting dimensions.
 */
void ALIGNMENT_Pushback(   ALIGNMENT* aln,
                           TRACE*     tr )
{
   /* if debugging, do edgechecks */
   #if DEBUG
   {
      /* if normal state, check bounds of normal dp matrix */
      /* if special state, check bounds of special dp matrix */
   }
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

/* compare two alignments */
int ALIGNMENT_Compare(  ALIGNMENT*     a,
                        ALIGNMENT*     b )
{
   return VECTOR_TRACE_Compare( a->traces, b->traces );
}

/*
 *  FUNCTION:  ALIGNMENT_Append()
 *  SYNOPSIS:  Append next state to Optimal Alignment.
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
inline
int ALIGNMENT_Append(   ALIGNMENT*   aln,       /* Traceback Alignment */
                        TRACE*       tr_ptr,    /* Traceback being Appended */
                        const int    st,        /* HMM state */
                        const int    q_0,       /* index in query/sequence */
                        const int    t_0 )      /* index in target/model */
{
   TRACE tr;

   /* for debugging: output traces as they are being added. */
   #if DEBUG 
   {
      // /* Add new state and (i,j) to trace */
      // int state_num[] = {MAT_ST, INS_ST, DEL_ST, SP_E, SP_N, SP_J, SP_C, SP_B, -1, -1};
      // if (st < 3) {
      //    fprintf( stderr, "%s:\t(%d,%d)\n", 
      //       STATE_NAMES[st], q_0, t_0 );
      // } 
      // else if (st >= 3 && st < 9) {
      //    fprintf( stderr, "%s:\t(%d,%d)\n", 
      //       STATE_NAMES[st], q_0, t_0 );
      // }
      // else {
      //    fprintf( stderr, "%s:\t(%d,%d)\n", 
      //       STATE_NAMES[st], q_0, t_0 );
      // }
   }
   #endif

   /* jump from current state to the prev state */
   switch (st)
   {
      /* Emit-on-Transition States: */
      case N_ST:
      case C_ST:
      case J_ST:
         tr.i = ( ( tr.st == st) ? q_0 : 0 );
         tr.j = 0;
         break;

      /* Non-Emitting States, not in Main Model: */
      case X_ST:
      case S_ST:
      case B_ST:
      case E_ST:
      case T_ST:
         tr.i = 0;
         tr.j = 0;
         break;

      /* Non-Emitting States, but in Main Model: */
      case D_ST:
         tr.i = 0;
         tr.j = t_0;
         break;

      /* Emitting States: */
      case M_ST:
      case I_ST:
         tr.i = q_0;
         tr.j = t_0;
         break;

      default:
         fprintf( stderr, "ERROR: Traceback failed. Invalid State Code occurred at [%d](%d,%d).\n", 
            st, q_0, t_0 );
         /* display valid states */
         for (int i = 0; i < NUM_ALL_STATES; i++) {
            fprintf( stderr, "[%d]:%s, ", i, STATE_NAMES[i] );
         }
         fprintf(stderr, "\n" );
         exit(EXIT_FAILURE);
   }

   tr.st = st;
   ALIGNMENT_Pushback( aln, &tr );

   return STATUS_SUCCESS;
}

/*
 *  FUNCTION:  ALIGNMENT_Find_Length()
 *  SYNOPSIS:  Scan <aln> for beginning, end, and length of alignments.
 */
void ALIGNMENT_Find_Length( ALIGNMENT* aln )
{
   /* scan traceback for all begin, end states */
   TRACE*   tr = aln->traces->data;
   int      N  = aln->traces->N;
   for (int i = 0; i < N; ++i) 
   {
      if ( tr[i].st == B_ST ) {
         VECTOR_INT_Pushback( aln->tr_beg, i + 1 );
      }
      if ( tr[i].st == E_ST ) {
         VECTOR_INT_Pushback( aln->tr_end, i - 1 );
      }
   }

   /* set beg,end vars to first found alignment (if at least one alignment) */
   if ( aln->tr_beg->N > 0 ) 
   {
      aln->beg       = aln->tr_beg->data[0];
      aln->end       = aln->tr_end->data[0];
      aln->aln_len   = aln->end - aln->beg + 1;
   }
}

/*
 *  FUNCTION:  ALIGNMENT_Reverse()
 *  SYNOPSIS:  Reverse order of <aln>.
 *             Alignments are built from back-to-front. This will correct that to normal ordering.
 */
void ALIGNMENT_Reverse( ALIGNMENT* aln )
{
   TRACE*   aln_0;   /* current trace */
   TRACE*   aln_1;   /* next trace */
   int      N;       /* alignment length */

   /* get length */
   N = aln->traces->N;

   /* update special states indexing */
   for (int i = 0; i < N; i++) 
   {
      aln_0 = &(aln->traces->data[i]);
      aln_1 = &(aln->traces->data[i+1]);

      if (aln_0->st == aln_1->st && (aln_0->st == N_ST || aln_0->st == C_ST || aln_0->st == J_ST) ) 
      {
         if (aln_0->i == 0 && aln_0->st > 0) 
         {
            aln_0->i = aln_1->i;
            aln_1->i = 0;
         }
      }
   }

   /* reverse state order */
   VECTOR_TRACE_Reverse( aln->traces );
}

/*
 *  FUNCTION:  ALIGNMENT_Build_HMMER_Style()
 *  SYNOPSIS:  Generate <aln> strings, HMMER-style.
 *             Expects <aln> has already been constructed.
 */
void ALIGNMENT_Build_HMMER_Style(   ALIGNMENT*     aln,
                                    SEQUENCE*      q_seq,
                                    HMM_PROFILE*   t_prof )
{
   /* alias for macros */
   HMM_PROFILE*   target      = t_prof;
   VECTOR_TRACE*  traceback   = aln->traces;
   TRACE*         tr          = NULL;

   int      i        = 0;
   int      beg      = aln->beg;
   int      end      = aln->end;
   int      aln_len  = end - beg + 1;

   /* allocate all string */
   aln->target_aln   = ERROR_realloc( aln->target_aln, sizeof(char) * (aln_len + 1) );
   aln->query_aln    = ERROR_realloc( aln->query_aln, sizeof(char) * (aln_len + 1) );
   aln->center_aln   = ERROR_realloc( aln->center_aln, sizeof(char) * (aln_len + 1) );
   aln->state_aln    = ERROR_realloc( aln->state_aln, sizeof(char) * (aln_len + 1) );

   /* counts */
   aln->num_matches  = 0;
   aln->num_misses   = 0;
   aln->num_gaps     = 0;

   /* construct consensus string if not already created */
   if ( t_prof->consensus == NULL ) {
      HMM_PROFILE_Set_Consensus( t_prof );
   }

   /* find beginnings and ends of traces */
   if ( aln->aln_len == -1 ) {
      ALIGNMENT_Find_Length( aln );
   }

   /* capture alignment (until END state) */
   int offset = 0;
   for (i = beg; i < end; i++, offset++)
   {
      tr          = &traceback->data[i];
      char t_ch   = t_prof->consensus[tr->j];
      char q_ch   = q_seq->seq[tr->i];
      char c_ch   = ' ';
      
      // printf("i=%d/%d; q_ch, t_ch, c_ch\n", i, end);

      /* output depends upon state */
      if ( tr->st == M_ST ) 
      { 
         aln->state_aln[offset]  = 'M';
         aln->target_aln[offset] = t_ch;
         aln->query_aln[offset]  = q_ch;
         /* if target and query match */
         if ( t_ch == q_ch ) {
            aln->num_matches++;
            aln->center_aln[offset] = t_ch;
         } 
         else {
            /* if target and query match (any case) */
            t_ch = tolower(t_ch);
            q_ch = tolower(q_ch);
            if ( t_ch == q_ch ) {
               aln->num_matches++;
               aln->center_aln[offset] = t_ch;
            }
            else {
               /* if score is positive */
               float msc;
               // msc = MSC(tr->j, q_ch);
               msc = 0;
               if ( MSC(tr->j, q_ch) > 0 ) {
                  aln->center_aln[offset] = '+';
               } 
               else {
                  aln->center_aln[offset] = '-';
               } 
            }
         }
      }
      if ( tr->st == I_ST ) 
      {
         aln->num_misses++;
         aln->state_aln[offset]  = 'I';
         /* insert corresponds to gap in target profile */
         aln->center_aln[offset] = ' ';
         aln->target_aln[offset] = '-';
         aln->query_aln[offset]  = q_ch;
      }
      if ( tr->st == D_ST ) {
         aln->num_misses++;
         aln->state_aln[offset]  = 'D';
         /* delete corresponds to gap in query sequence */
         aln->center_aln[offset] = ' ';
         aln->target_aln[offset] = t_ch;
         aln->query_aln[offset]  = '-';
      }
   }

   aln->state_aln[offset]     = '\0';
   aln->center_aln[offset]    = '\0';
   aln->target_aln[offset]    = '\0';
   aln->query_aln[offset]     = '\0';

   aln->perc_id = (float)aln->num_matches / (float)aln->traces->N;

   // printf("state:  %s\n", aln->state_aln );
   // printf("target: %s\n", aln->target_aln );
   // printf("center: %s\n", aln->center_aln );
   // printf("query:  %s\n", aln->query_aln );
}

/*
 *  FUNCTION:  ALIGNMENT_Build_MMSEQS_Style()
 *  SYNOPSIS:  Generate <aln> strings, HMMER-style (Cigar-style).
 *             Expects <aln> has already been constructed.
 *   EXAMPLE:
 */
void ALIGNMENT_Build_MMSEQS_Style(  ALIGNMENT*     aln,
                                    SEQUENCE*      query,
                                    HMM_PROFILE*   target )
{
   VECTOR_TRACE*  traceback   = aln->traces;

   /* cigar alignment can't be longer than twice the size of the alignment */
   aln->cigar_aln = ERROR_malloc( sizeof(char) * ((2 * aln->traces->N) + 1) );
   aln->cigar_aln[0] = '\0';
   /* place to cast run length and state to string */
   char     cigar_buffer[32];
   /* offset into cigar */
   int      offset      = 0;
   /* tracks number of consequetive states */
   int      run_len     = 0;

   TRACE*   tr          = NULL;
   int      prv_st      = -1;
   int      cur_st      = -1;

   /* unrolled first iteration */
   tr       = &traceback->data[0];
   prv_st   = tr->st;

   for (int i = 1; i < traceback->N; i++)
   {
      tr       = &traceback->data[i];
      cur_st   = tr->st;
      if ( cur_st == prv_st ) {
         run_len += 1;
      } else {
         sprintf( cigar_buffer, "%s%d", STATE_CHARS[prv_st], run_len );
         strcat( aln->cigar_aln, cigar_buffer );
         run_len = 1;
      }
      prv_st = cur_st;
   }

   printf("Cigar: %s\n", aln->cigar_aln);
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
