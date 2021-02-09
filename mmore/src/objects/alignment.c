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
#include "../utilities/_utilities.h"
#include "_objects.h"

/* header */
#include "alignment.h"

/*! FUNCTION:  ALIGNMENT_Create()
 *  SYNOPSIS:  Constructs <aln>, allocates memory and returns pointer.
 */
ALIGNMENT* 
ALIGNMENT_Create()
{
   ALIGNMENT *aln       = NULL;
   const int min_size   = 256;

   aln = ERROR_malloc( sizeof(ALIGNMENT) );

   aln->seq_beg      = VECTOR_INT_Create();
   aln->seq_end      = VECTOR_INT_Create();
   aln->sequence     = VECTOR_CHAR_Create();

   aln->tr_beg       = VECTOR_INT_Create();
   aln->tr_end       = VECTOR_INT_Create();
   aln->traces       = VECTOR_TRACE_Create();

   aln->Q            = -1;
   aln->T            = -1;
   aln->aln_len      = 0;

   aln->perc_id      = 0.0;
   aln->num_matches  = 0;
   aln->num_gaps     = 0;
   aln->num_misses   = 0;

   aln->cigar_aln    = VECTOR_CHAR_Create();
   aln->is_cigar_aln = false;
   aln->target_aln   = VECTOR_CHAR_Create();
   aln->center_aln   = VECTOR_CHAR_Create();
   aln->query_aln    = VECTOR_CHAR_Create();
   aln->state_aln    = VECTOR_CHAR_Create();
   aln->is_hmmer_aln = false;

   ALIGNMENT_Resize(aln, min_size);

   return aln;
}

/*! FUNCTION:  ALIGNMENT_Destroy()
 *  SYNOPSIS:  Destroys <aln>, frees memory and return NULL pointer.
 */
ALIGNMENT* 
ALIGNMENT_Destroy(   ALIGNMENT* aln )
{
   if (aln == NULL) return aln;
   
   aln->seq_beg      = VECTOR_INT_Destroy( aln->seq_beg );
   aln->seq_end      = VECTOR_INT_Destroy( aln->seq_end );
   aln->sequence     = VECTOR_CHAR_Destroy( aln->sequence );

   aln->tr_beg       = VECTOR_INT_Destroy( aln->tr_beg );
   aln->tr_end       = VECTOR_INT_Destroy( aln->tr_end );
   aln->traces       = VECTOR_TRACE_Destroy( aln->traces );

   aln->cigar_aln    = VECTOR_CHAR_Destroy( aln->cigar_aln );
   aln->target_aln   = VECTOR_CHAR_Destroy( aln->target_aln );
   aln->center_aln   = VECTOR_CHAR_Destroy( aln->center_aln );
   aln->query_aln    = VECTOR_CHAR_Destroy( aln->query_aln );
   aln->state_aln    = VECTOR_CHAR_Destroy( aln->state_aln );

   aln = ERROR_free( aln );

   return aln;
}

/*! FUNCTION:  ALIGNMENT_Reuse()
 *  SYNOPSIS:  Wipes <aln>'s old data for reuse, sets dimensions <Q x T>.
 */
void 
ALIGNMENT_Reuse(  ALIGNMENT*  aln,
                  int         Q,
                  int         T )
{
   VECTOR_INT_Reuse( aln->seq_beg );
   VECTOR_INT_Reuse( aln->seq_end );
   VECTOR_CHAR_Reuse( aln->sequence );

   VECTOR_INT_Reuse( aln->tr_beg );
   VECTOR_INT_Reuse( aln->tr_end );
   VECTOR_TRACE_Reuse( aln->traces );

   VECTOR_CHAR_Reuse( aln->cigar_aln );
   aln->is_cigar_aln = false;
   VECTOR_CHAR_Reuse( aln->target_aln );
   VECTOR_CHAR_Reuse( aln->center_aln );
   VECTOR_CHAR_Reuse( aln->query_aln );
   VECTOR_CHAR_Reuse( aln->state_aln );
   aln->is_hmmer_aln = false;
}

/*! FUNCTION:  RESULTS_Reuse()
 *  SYNOPSIS:  Push trace onto end of alignment.
 */
void 
ALIGNMENT_Pushback(     ALIGNMENT*  aln,
                        TRACE*      tr )
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

/*! FUNCTION:  ALIGNMENT_Resize()
 *  SYNOPSIS:  Resize <aln>'s trace array to <size>.
 */
void 
ALIGNMENT_Resize(    ALIGNMENT*   aln,
                     size_t       size )
{
   VECTOR_TRACE_Resize( aln->traces, size );
}

/*! FUNCTION:  ALIGNMENT_GrowTo()
 *  SYNOPSIS:  Resize <aln>'s trace array to <size> only if array is smaller than current size.
 */
void 
ALIGNMENT_GrowTo(    ALIGNMENT*   aln,
                     size_t       size )
{
   VECTOR_TRACE_GrowTo( aln->traces, size );
}

/*! FUNCTION:  ALIGNMENT_Resize()
 *  SYNOPSIS:  Compare alignments <a> and <b>.
 *  RETURN:    Zero if equal.
 */
int 
ALIGNMENT_Compare(   ALIGNMENT*     a,
                     ALIGNMENT*     b )
{
   return VECTOR_TRACE_Compare( a->traces, b->traces );
}

/*! FUNCTION:  ALIGNMENT_Append()
 *  SYNOPSIS:  Append trace <st_cur, q_0, t_0> to <aln>.
 */
inline
int 
ALIGNMENT_Append(    ALIGNMENT*   aln,       /* Traceback Alignment */
                     const int    st,        /* HMM state */
                     const int    q_0,       /* index in query/sequence */
                     const int    t_0 )      /* index in target/model */
{
   TRACE    tr, prv_tr;

   /* jump from current state to the prev state */
   switch (st)
   {
      /* Emit-on-Transition States: */
      case N_ST:
      case C_ST:
      case J_ST:
         prv_tr = VEC_X( aln->traces, aln->traces->N - 1 );
         tr.q_0 = ( ( tr.st == prv_tr.st ) ? q_0 : 0 );
         tr.t_0 = 0;
         break;

      /* Non-Emitting States, not in Main Model: */
      case X_ST:
      case S_ST:
      case B_ST:
      case E_ST:
      case T_ST:
         tr.q_0 = 0;
         tr.t_0 = 0;
         break;

      /* Non-Emitting States, but in Main Model: */
      case D_ST:
         tr.q_0 = 0;
         tr.t_0 = t_0;
         break;

      /* Emitting States: */
      case M_ST:
      case I_ST:
         tr.q_0 = q_0;
         tr.t_0 = t_0;
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

/*! FUNCTION:  ALIGNMENT_Find_Length()
 *  SYNOPSIS:  Scan <aln> for beginning, end, and length of alignment from <B> to <E> states. 
 *             Stores <beg> and <end>.
 */
void 
ALIGNMENT_Find_Length(  ALIGNMENT*  aln )
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

   /* set beg, end vars to first found alignment (if at least one alignment) */
   if ( aln->tr_beg->N > 0 ) 
   {
      aln->beg       = aln->tr_beg->data[0];
      aln->end       = aln->tr_end->data[0];
      aln->aln_len   = aln->end - aln->beg + 1;
   }
}

/*! FUNCTION:  ALIGNMENT_SetEndpoints()
 *  SYNOPSIS:  Sets <beg> and <end> endpoint indexes of the <aln> alignment.
 */
void 
ALIGNMENT_SetEndpoints(    ALIGNMENT*  aln,
                           int         beg,
                           int         end )
{
   aln->beg       = beg;
   aln->end       = end;
   aln->aln_len   = end - beg;
}

/*! FUNCTION:  ALIGNMENT_Reverse()
 *  SYNOPSIS:  Reverse order of <aln>.
 *             For use with alignments built from back-to-front via backtrace. 
 *             This will correct that to normal ordering.
 */
void 
ALIGNMENT_Reverse(   ALIGNMENT*  aln )
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
         if (aln_0->q_0 == 0 && aln_0->st > 0) 
         {
            aln_0->q_0 = aln_1->q_0;
            aln_1->q_0 = 0;
         }
      }
   }

   /* reverse state order */
   VECTOR_TRACE_Reverse( aln->traces );
}

/*! FUNCTION:  ALIGNMENT_Build_HMMER_Style()
 *  SYNOPSIS:  Generate <aln> strings, HMMER-style. 
 *             Stores in <target_aln>, <center_aln>, <query_aln>.
 *             Expects <aln> has already been constructed.
 */
void 
ALIGNMENT_Build_HMMER_Style(  ALIGNMENT*     aln,
                              SEQUENCE*      q_seq,
                              HMM_PROFILE*   t_prof )
{
   /* alias for macros */
   HMM_PROFILE*   target      = t_prof;   /* need alias for MSC() to work */
   VECTOR_TRACE*  traceback   = aln->traces;
   TRACE*         tr          = NULL;
   /* query and target consensus sequences as strings */
   STR            qseq        = SEQUENCE_GetSeq( q_seq );
   STR            tseq        = HMM_PROFILE_GetConsensus( t_prof );
   int            T           = t_prof->N;
   int            Q           = q_seq->N;
   /* loop vars */
   int            i           = 0;
   int            i_beg       = aln->beg;
   int            i_end       = aln->end;
   int            aln_len     = i_end - i_beg + 1;
   char           t_ch;
   char           q_ch;
   char           c_ch;
   float          msc;
   /* special characters for alignment */
   const char     pos_ch   = '+';
   const char     neg_ch   = '-';
   const char     gap_ch   = '-';
   
   /* counters */
   aln->num_matches  = 0;
   aln->num_misses   = 0;
   aln->num_gaps     = 0;
   /* streak counter */
   int match_streak  = 0;

   /* allocate for entire string */
   VECTOR_CHAR_Reuse( aln->target_aln );
   VECTOR_CHAR_Reuse( aln->query_aln );
   VECTOR_CHAR_Reuse( aln->center_aln );
   VECTOR_CHAR_Reuse( aln->state_aln );
   VECTOR_CHAR_GrowTo( aln->target_aln, aln_len + 1 );
   VECTOR_CHAR_GrowTo( aln->query_aln, aln_len + 1 );
   VECTOR_CHAR_GrowTo( aln->center_aln, aln_len + 1 );
   VECTOR_CHAR_GrowTo( aln->state_aln, aln_len + 1 );

   /* find beginnings and ends of traces */
   if ( aln->aln_len == -1 ) {
      ALIGNMENT_Find_Length( aln );
   }

   /* capture alignment (until END state) */
   int offset = 0;
   for (i = i_beg; i <= i_end; i++, offset++)
   {
      /* get emitted residue at position in the alignment */
      tr       = &VEC_X( traceback, i );
      t_ch     = STR_GetChar( tseq, tr->t_0 );
      q_ch     = STR_GetChar( qseq, tr->q_0 );
      c_ch     = ' ';

      /* center symbol depends upon state */
      if ( tr->st == M_ST ) 
      { 
         match_streak++;
         VECTOR_CHAR_Set( aln->state_aln,    offset, 'M' );
         VECTOR_CHAR_Set( aln->target_aln,   offset, t_ch );
         VECTOR_CHAR_Set( aln->query_aln,    offset, q_ch );

         /* center alignment depends on match score between target and query: */
         /* if target and query match (same case) */
         if ( t_ch == q_ch ) {
            aln->num_matches++;
            VECTOR_CHAR_Set( aln->center_aln, offset, t_ch );
         } 
         /* else, if target and query match (any case) */
         else {
            t_ch = tolower( t_ch );
            q_ch = tolower( q_ch );
            if ( t_ch == q_ch ) {
               aln->num_matches++;
               VECTOR_CHAR_Set( aln->center_aln, offset, t_ch );
            }
            /* else, if score is positive or negative */
            else {
               /* get match score at current position in alignment */
               msc = MSC( tr->t_0, AA_REV[q_ch] );
               /* if match score is positive */
               if ( msc > 0 ) {
                  VECTOR_CHAR_Set( aln->center_aln, offset, '+' );
               } 
               /* if match score is negative */
               else {
                  VECTOR_CHAR_Set( aln->center_aln, offset, '-' );
               } 
            }
         }
      }
      elif ( tr->st == I_ST ) 
      {
         /* if we are leaving a matching region, add a gap */
         if (match_streak > 0) {
            aln->num_gaps++;
         }
         match_streak = 0;
         aln->num_misses++;
         /* insert corresponds to gap in target profile */
         VECTOR_CHAR_Set( aln->state_aln,    offset, 'I' );
         VECTOR_CHAR_Set( aln->center_aln,   offset, ' ' );
         VECTOR_CHAR_Set( aln->target_aln,   offset, '-' );
         VECTOR_CHAR_Set( aln->query_aln,    offset, q_ch );
      }
      elif ( tr->st == D_ST ) 
      {
         /* if we are leaving a matching region, add a gap */
         if (match_streak > 0) {
            aln->num_gaps++;
         }
         match_streak = 0;
         aln->num_misses++;
         /* delete corresponds to gap in query sequence */
         VECTOR_CHAR_Set( aln->state_aln,    offset, 'D' );
         VECTOR_CHAR_Set( aln->center_aln,   offset, ' ' );
         VECTOR_CHAR_Set( aln->target_aln,   offset, t_ch );
         VECTOR_CHAR_Set( aln->query_aln,    offset, '-' );
      }
   }
   /* terminate strings with null termination character */
   VECTOR_CHAR_Set( aln->state_aln,    offset, NULL_CHAR );
   VECTOR_CHAR_Set( aln->center_aln,   offset, NULL_CHAR );
   VECTOR_CHAR_Set( aln->target_aln,   offset, NULL_CHAR );
   VECTOR_CHAR_Set( aln->query_aln,    offset, NULL_CHAR );

   /* percent identity is percentage of alignment that is a direct match */
   aln->perc_id = (float)aln->num_matches / (float)aln->traces->N;

   aln->is_hmmer_aln = true;
}

/*! FUNCTION:  ALIGNMENT_Build_MMSEQS_Style()
 *  SYNOPSIS:  Generate <aln> strings, MMSEQS-style. 
 *             Stores in <cigar_aln>.
 *             Expects <aln> has already been constructed.
 */
void 
ALIGNMENT_Build_MMSEQS_Style(    ALIGNMENT*     aln,
                                 SEQUENCE*      query,
                                 HMM_PROFILE*   target )
{
   VECTOR_TRACE*  traceback   = aln->traces;

   /* cigar alignment can't be longer than twice the size of the alignment */
   VECTOR_CHAR_Reuse( aln->cigar_aln );

   /* place to cast run length and state to string */
   char     cigar_buffer[64];    /* buffer for building string */
   int      offset      = 0;     /* offset into cigar string */
   int      run_len     = 0;     /* tracks number of consequetive states */

   TRACE*   tr          = NULL;  /* current trace */
   int      prv_st      = -1;    /* previous state */
   int      cur_st      = -1;    /* current state */

   /* unrolled first iteration */
   tr       = &VEC_X( traceback, 0 );
   prv_st   = tr->st;
   run_len  = 1;

   for (int i = 1; i < traceback->N; i++)
   {
      tr       = &VEC_X( traceback, i );
      cur_st   = tr->st;

      /* if cuurent state is the same as previous, continue run */
      if ( cur_st == prv_st ) {
         run_len += 1;
      }
      /* if cuurent state change from previous, add to alignment */ 
      else {
         /* format cigar string addition to {ST_CHAR}{RUN_LENGTH} e.g. M12, D4, B1 */
         sprintf( cigar_buffer, "%s%d", STATE_CHARS[prv_st], run_len );
         VECTOR_CHAR_Append( aln->cigar_aln, cigar_buffer, strlen(cigar_buffer) );
         run_len = 1;
      }
      prv_st = cur_st;
   }
   /* add last run to alignment */
   sprintf( cigar_buffer, "%s%d", STATE_CHARS[prv_st], run_len );
   VECTOR_CHAR_Append( aln->cigar_aln, cigar_buffer, strlen(cigar_buffer) );

   /* add null character to end of string */
   VECTOR_CHAR_Pushback( aln->cigar_aln, NULL_CHAR );

   aln->is_cigar_aln = true;
}

/*! FUNCTION:  ALIGNMENT_Dump()
 *  SYNOPSIS:  Outputs <aln> to open file pointer <fp>.
 */
void 
ALIGNMENT_Dump(   ALIGNMENT*  aln,
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
      fprintf(fp, "[%d](%s,%d,%d)\n", i, STATE_NAMES[tr[i].st], tr[i].q_0, tr[i].t_0);
   }
}

/*! FUNCTION:  ALIGNMENT_Save()
 *  SYNOPSIS:  Save <aln> to file at location <filename>. 
 *             Handles opening and closing of file.
 */
void 
ALIGNMENT_Save(   ALIGNMENT*  aln,
                  char*       filename )
{
   FILE* fp = fopen(filename, "w");
   ALIGNMENT_Dump(aln, fp);
   fclose(fp);
}


