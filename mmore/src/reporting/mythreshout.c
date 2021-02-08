/*******************************************************************************
 *  FILE:      myout.c
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
#include "myout.h"

/* === MYOUT FUNCTIONS === */

/*    FUNCTION:   REPORT_mythreshout_header()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after m8, see example)
 */
void REPORT_mythreshout_header(  WORKER*  worker,
                                 FILE*    fp )
{
   const int num_fields = 8;
   const char* headers[] = {
      "target-hmm",
      "query-seq",
      "viterbi",
      "cloud-max",
      "cloud-compo",
      "bound-max",
      "dom-max",
      "dom-compo"
   };

   REPORT_header(fp, headers, num_fields);
}

/*    FUNCTION:   REPORT_mythreshout_entry()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
void REPORT_mythreshout_entry(   WORKER*  worker,
                                 RESULT*  result,
                                 FILE*    fp )
{
   HMM_PROFILE*   t_prof   = worker->t_prof;
   SEQUENCE*      q_seq    = worker->q_seq;
   ALL_SCORES*    scores   = &result->scores;
   SCORES*        final    = &result->final_scores;

   const int num_fields = 8;
   const int sig_digits = 3;

   const GEN fields[] = {
      GEN_Wrap( &t_prof->name,                    DATATYPE_STRING,  sizeof(char*) ),
      GEN_Wrap( &q_seq->name,                     DATATYPE_STRING,  sizeof(char*) ),
      GEN_Wrap( &final->viterbi_eval,             DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Wrap( &scores->threshold_cloud_max,     DATATYPE_FLOAT,   sizeof(float) ),  
      GEN_Wrap( &scores->threshold_cloud_compo,   DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Wrap( &scores->threshold_bound_max,     DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Wrap( &scores->threshold_dom_max,       DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Wrap( &scores->threshold_dom_compo,     DATATYPE_FLOAT,   sizeof(float) )
   };

   REPORT_entry( fp, fields, num_fields, sig_digits );
}

/*    FUNCTION:   REPORT_mythreshout_footer()
 *    SYNOPSIS:   Print footer
 *                (modeled after HMMER, see example)              
 */
void REPORT_mythreshout_footer(    WORKER*  worker,
                                    FILE*    fp )
{
   fprintf( fp, "# [ok] [myout]\n" );
}