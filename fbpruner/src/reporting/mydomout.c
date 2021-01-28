/*******************************************************************************
 *  FILE:      mydomout.c
 *  PURPOSE:   Reporting Subroutines for generating mydomout.
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
#include "mydomout.h"

/* === MYDOMOUT FUNCTIONS === */

/* === MYDOM OUTPUT === */
/* Custom-style output for domain-specific scoring 
   Description: The file is formatted as a tab-separated list with 12 columns: 
   - (1) query profile identifier
   - (2) target sequences identifier 
   - (3) e-value, 
   - (4) bitscore before composition bias correction
   - (5) composition bias score
   - (6) bitscore after composition bias correction
   - (7) mmseqs viterbi bitscore
   - (8) total DP matrix cells used by full viterbi/forward-backward
   - (9) DP matrix cloud of cells used by MMORE
   - (10) start-end range of query in cloud
   - (11) start-end range of target in cloud
   - (12) time to run given search
 */

/* EXAMPLE:
 *
   [header/]
   #  qid          tid                                 %id    aln   mis  gap  qb   qe   tb    te    eval       bit 
   #  ---          ---                                 ---    ---   ---  ---  ---  ---  ---   ---   ----       ---
   [entry/]
      Aldolase_II  Aldolase_II/544/29-203/211-382      0.252  151   105  0    21   171  51    191   2.631E-07  54
      Aldolase_II  decoy94676                          0.230  91    69   0    66   156  164   253   2.572E+01  29
      Aldolase_II  decoy95461                          0.553  24    10   0    16   39   102   124   4.583E+01  29
   [footer/]
   # [ok]
 *
 */

/*    FUNCTION:   REPORT_mydomout_header()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
void REPORT_mydomout_header(  WORKER*  worker,
                              FILE*    fp )
{
   fprintf( fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
      "query",
      "target",
      "n_domains",
      "evalue",
      "pre-sc",
      "comp-bias",
      "seq-sc",
      "dom-sum-sc",
      "mmseqs-vit-sc",
      "total-cells",
      "MMORE-cells",
      "perc-cells",
      "q-range",
      "t-range",
      "time"
   );
   fprintf( fp, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
      "------",
      "------",
      "---------",
      "------",
      "------",
      "--------",
      "------",
      "----------",
      "-------------",
      "----------",
      "----------",
      "----------",
      "-------",
      "-------",
      "----"
   );
}

/*    FUNCTION:   REPORT_mydomout_entry()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
void REPORT_mydomout_entry(   WORKER*  worker,
                              RESULT*  result,
                              FILE*    fp )
{
   HMM_PROFILE*   t_prof   = worker->t_prof;
   SEQUENCE*      q_seq    = worker->q_seq;
   ALIGNMENT*     aln      = worker->traceback;
   DOMAIN_DEF*    dom_def  = worker->dom_def;

   /* if no domains were found, skip entry */
   if ( dom_def->n_domains <= 0 ) {
      return;
   }

   int      best_idx    = dom_def->best;
   RANGE    dom_rng     = dom_def->best_range;
   float    dom_fwdsc   = dom_def->best_fwdsc;
   float    dom_presc   = dom_def->best_presc;
   float    dom_bias    = dom_def->best_bias;
   float    dom_sc      = dom_def->best_sc;
   float    dom_sumsc   = dom_def->dom_sumsc;

   float percent_cells = (float)result->cloud_cells / (float)result->total_cells;

   fprintf( fp, "%s\t%s\t%d/%d\t%.3g\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%.5f\t%d-%d\t%d-%d\t%.5f\n", 
      t_prof->name,                    /* target name */
      q_seq->name,                     /* query name */
      best_idx+1, dom_def->n_domains,  /* domain id */
      result->final_scores.eval,       /* evalue */
      dom_presc,                       /* seq scores before correction */
      dom_bias,                        /* seq bias */
      dom_sc,                          /* seq score after correction */
      dom_sumsc,                       /* sum of all domain scores */
      result->vit_natsc,               /* viterbi score (in mmore, this comes from mmseqs) */
      result->total_cells,             /* total number of cells computed by full viterbi */
      result->cloud_cells,             /* total number of cells computed by mmore */
      percent_cells,                   /* percent of total cells computed by mmore */
      result->target_start,            /* start of target range */
      result->target_end,              /* end of target range */
      dom_rng.beg,                     /* start of query range */
      dom_rng.end,                     /* end of query range */
      result->time                     /* time for entire iteration */
   );
}

/*    FUNCTION:   REPORT_mydomout_footer()
 *    SYNOPSIS:   Print footer
 *                (modeled after HMMER, see example)              
 */
void REPORT_mydomout_footer(    WORKER*  worker,
                              FILE*    fp )
{
   fprintf( fp, "# [ok] [myout]\n" );
}

