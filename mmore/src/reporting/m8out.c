/*******************************************************************************
 *  FILE:      report.c
 *  PURPOSE:   Reporting Subroutines for generating output.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       - 
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
#include "../objects/structs.h"
#include "../utilities/_utilities.h"
#include "../objects/_objects.h"
#include "../parsers/_parsers.h"
#include "../algs_linear/_algs_linear.h"
#include "../algs_quad/_algs_quad.h"
#include "../algs_naive/_algs_naive.h"

/* header */
#include "_reporting.h"
#include "m8out.h"

/* === M8 OUTPUT === */
/* BLAST-style output 
   Description: The file is formatted as a tab-separated list with 12 columns: 
   - (1,2) identifiers query and target sequences/profiles, 
   - (3) sequence identity, 
   - (4) alignment leng
   - (5) number of mismatches, 
   - (6) number of gap openings, 
   - (7,8) domain (start,end)-position in query and in target, 
   - (9,10) domain
   - (11) E-value, and 
   - (12) bit score.
 */

/* EXAMPLE:
 *
   [header/]
   #  qid         tid                                  %id    aln   mis  gap  qb   qe   tb    te    eval       bit 
   #  ----------  -------------------------------      ---    ---   ---  ---  ---  ---  ---   ---   ----       ---
   [entry/]
      Aldolase_II  Aldolase_II/544/29-203/211-382      0.252  151   105  0    21   171  51    191   2.631E-07  54
      Aldolase_II  decoy94676                          0.230  91    69   0    66   156  164   253   2.572E+01  29
      Aldolase_II  decoy95461                          0.553  24    10   0    16   39   102   124   4.583E+01  29
   [footer/]
   # [ok]
 *
 */

/*!   FUNCTION:   REPORT_m8out_header()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
STATUS_FLAG 
REPORT_m8out_header(    WORKER*  worker,
                        FILE*    fp )
{
   const int num_fields = 13;
   const char* headers[] = {
      "query-hmm",
      "target-seq",
      "pident",
      "alnlen",
      "mismatch",
      "gapopen",
      "qstart",
      "qend",
      "tstart",
      "tend",
      "eval",
      "bits",
      "cigar"
   };

   REPORT_header(fp, headers, num_fields);
}

/*!   FUNCTION:   REPORT_m8out_entry()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
STATUS_FLAG 
REPORT_m8out_entry(     WORKER*  worker,
                        RESULT*  result,
                        FILE*    fp )
{
   HMM_PROFILE*   t_prof      = worker->t_prof;
   SEQUENCE*      q_seq       = worker->q_seq;
   ALIGNMENT*     aln         = worker->trace_vit;
   STR            cigar_aln   = NULL;

   TRACE*         beg         = &VEC_X( aln->traces, aln->beg );
   TRACE*         end         = &VEC_X( aln->traces, aln->end );

   if ( aln->is_cigar_aln == false ) {
      ALIGNMENT_Build_MMSEQS_Style( worker->trace_vit, worker->q_seq, worker->t_prof );
   }
   cigar_aln = VECTOR_CHAR_GetArray( aln->cigar_aln );
   cigar_aln = ( STR_GetLength(cigar_aln) > 0 ? cigar_aln : "--" );

   fprintf( fp, "%s\t%s\t%.3f\t%ld\t%d\t%d\t%d\t%d\t%d\t%d\t%9.2g\t%6.1f\t%s\n", 
      t_prof->name,                    /* target name */
      q_seq->name,                     /* query name */
      aln->perc_id,                    /* percent id (number matches) */
      aln->aln_len,                    /* alignment length */
      aln->num_misses,                 /* number of mismatches */    
      aln->num_gaps,                   /* number of gap openings */
      beg->t_0,                        /* query start */
      end->t_0,                        /* query end */
      beg->q_0,                        /* target start */
      end->q_0,                        /* target end */
      result->final_scores.eval,       /* evalue */
      result->final_scores.nat_sc,     /* bits */
      cigar_aln                        /* MMSEQS-style, cigar alignment */
   );
}

/*!   FUNCTION:   REPORT_m8out_footer()
 *    SYNOPSIS:   Print footer
 *                (modeled after HMMER, see example)              
 */
STATUS_FLAG 
REPORT_m8out_footer(    WORKER*  worker,
                        FILE*    fp )
{
   fprintf( fp, "# [ok] [m8out]\n" );
}
