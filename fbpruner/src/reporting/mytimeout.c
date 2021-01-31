/*******************************************************************************
 *  FILE:      myout.c
 *  PURPOSE:   Reporting functions for generating mytimeout format output.
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
#include "mytimeout.h"

/* === MYTIMEOUT FUNCTIONS === */

/* === MYTIMEOUT OUTPUT === */
/* Custom-style output for logging times for subroutines in pipeline.
   Description: The file is formatted as a tab-separated list with the following columns: 
   - (1) total time
   - (2) target (hmm model) load time (zero if it doesn't need to load new model)
   - (3) query (sequence) load time
   - (4) cloud search fwd time
   - (5) cloud search bck time
   - (6) cloud union time
   - (7) cloud reorient time
   - (8) pruned fwd linear time
   - (9) pruned bck linear time
   - (10) sparse matrix build time
   - (11) pruned fwd sparse time
   - (12) pruned bck sparse time
   - (13) posterior sparse time
   - (14) bias-corr sparse time
   - (15) optimal alignment sparse time
 */

/*!  FUNCTION:    REPORT_mytimeout_header()
 *   SYNOPSIS:    Print header for all mytimeout fields.
 */
STATUS_FLAG 
REPORT_mytimeout_header(   WORKER*  worker,
                           FILE*    fp )
{
   const int   num_fields = 23;
   const char* headers[] = {
      "target-hmm",
      "query-seq",
      "total",
      "hmm-load",
      "seq-load",
      "cloud-fwd",
      "cloud-bck",
      "cloud-union",
      "cloud-reorient",
      "bound-fwd-lin",
      "bound-bck-lin",
      "sparse-build",
      "bound-fwd-sp",
      "bound-bck-sp",
      "posterior-sp",
      "decodedom-sp",
      "biascorr-sp",
      "optacc-sp",
      "bound-fwd-dom",
      "bound-bck-dom",
      "posterior-dom",
      "biascorr-dom",
      "optacc-dom"
   };

   REPORT_header( fp, headers, num_fields );

   return STATUS_SUCCESS;
}

/*!  FUNCTION:    REPORT_myout_entry()
 *   SYNOPSIS:    Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
STATUS_FLAG 
REPORT_mytimeout_entry(    WORKER*  worker,
                           RESULT*  result,
                           FILE*    fp )
{
   HMM_PROFILE*   t_prof         = worker->t_prof;
   SEQUENCE*      q_seq          = worker->q_seq;
   TIMES*         times          = worker->times;
   TIMES*         time_totals    = worker->times_totals;

   const int num_fields = 23;
   const int sig_digits = 7;

   const GEN fields[] = {
      GEN_Create( &t_prof->name,         DATATYPE_STRING,  sizeof(char*) ),
      GEN_Create( &q_seq->name,          DATATYPE_STRING,  sizeof(char*) ),
      GEN_Create( &times->total,         DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Create( &times->load_target,   DATATYPE_FLOAT,   sizeof(float) ),  
      GEN_Create( &times->load_query,    DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Create( &times->lin_cloud_fwd, DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Create( &times->lin_cloud_bck, DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Create( &times->lin_merge,     DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Create( &times->lin_reorient,  DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Create( &times->lin_bound_fwd, DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Create( &times->lin_bound_bck, DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Create( &times->sp_build_mx,   DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Create( &times->sp_bound_fwd,  DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Create( &times->sp_bound_bck,  DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Create( &times->sp_posterior,  DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Create( &times->sp_decodedom,  DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Create( &times->sp_biascorr,   DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Create( &times->sp_optacc,     DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Create( &times->dom_bound_fwd, DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Create( &times->dom_bound_bck, DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Create( &times->dom_posterior, DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Create( &times->dom_biascorr,  DATATYPE_FLOAT,   sizeof(float) ),
      GEN_Create( &times->dom_optacc,    DATATYPE_FLOAT,   sizeof(float) ),
   };

   REPORT_entry( fp, fields, num_fields, sig_digits );
}

/*!  FUNCTION:    REPORT_myout_footer()
 *   SYNOPSIS:    Print footer.  Gives [ok] to verify file program executed succesfully. And file format type.          
 */
STATUS_FLAG 
REPORT_mytimeout_footer(   WORKER*  worker,
                           FILE*    fp )
{
   fprintf( fp, "# [ok] [mytimeout]\n" );
}

/*!  FUNCTION:    REPORT_mytimeout_totals()
 *   SYNOPSIS:    Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
STATUS_FLAG 
REPORT_mytimeout_totals(   WORKER*  worker,
                           RESULT*  result,
                           FILE*    fp )
{
   HMM_PROFILE*   t_prof   = worker->t_prof;
   SEQUENCE*      q_seq    = worker->q_seq;
   TIMES*         times    = worker->times_totals;

   const int   num_fields  = 23;
   const int   sig_digits  = 7;
   const int   perc_digits = 5;
   const int   pad         = 24;

   char* headers[] = {
      "program",
      "program-outside-loop",
      "loop-total",
      "loop-sum",
      "loop-unaccounted",
      "hmm-load",
      "seq-load",
      "cloud-fwd",
      "cloud-bck",
      "cloud-union",
      "cloud-reorient",
      "bound-fwd-lin",
      "bound-bck-lin",
      "sparse-build",
      "bound-fwd-sp",
      "bound-bck-sp",
      "posterior-sp",
      "decodedom-sp",
      "biascorr-sp",
      "optacc-sp",
      "bound-fwd-dom",
      "bound-bck-dom",
      "posterior-dom",
      "biascorr-dom",
      "optacc-dom"
   };

   float data[] = {
      times->program_runtime,
      times->program_runtime - times->total,
      times->total,
      0.0f,
      times->total,
      times->load_target,   
      times->load_query,  
      times->lin_cloud_fwd,  
      times->lin_cloud_bck,  
      times->lin_merge,     
      times->lin_reorient,  
      times->lin_bound_fwd,  
      times->lin_bound_bck, 
      times->sp_build_mx,  
      times->sp_bound_fwd, 
      times->sp_bound_bck, 
      times->sp_posterior,
      times->sp_decodedom, 
      times->sp_biascorr,  
      times->sp_optacc,
      times->dom_bound_fwd,
      times->dom_bound_bck,
      times->dom_posterior,
      times->dom_biascorr,
      times->dom_optacc
   };

   /* sum subroutines */
   for (int i = 5; i < 25; i++) {
      data[3] += data[i];
      data[4] -= data[i];
   }

   /* headers */
   fprintf(fp, "#%*s:     %*s    %*s\n", 
      pad-1, "TOTALS", 
      sig_digits+5, "TIME", 
      perc_digits+5, "PERCENT" );
   fprintf(fp, "#%*s:     %*s    %*s\n", 
      pad-1, "------", 
      sig_digits+5, "----", 
      perc_digits+5, "-------" );

   /* data rows */
   for (int i = 0; i < 2; i++) {
      fprintf(fp, "%*s:     %*.*f    %*.*f%%\n", 
         pad, headers[i], 
         sig_digits+5, sig_digits, data[i], 
         perc_digits+5, perc_digits, data[i]/data[0] );
   }
   fprintf(fp, "\n");
   for (int i = 2; i < 5; i++) {
      fprintf(fp, "%*s:     %*.*f    %*.*f%%\n", 
         pad, headers[i], 
         sig_digits+5, sig_digits, data[i], 
         perc_digits+5, perc_digits, data[i]/data[0] );
   }
   fprintf(fp, "\n");
   for (int i = 5; i < 25; i++) {
      fprintf(fp, "%*s:     %*.*f    %*.*f%%\n", 
         pad, headers[i], 
         sig_digits+5, sig_digits, data[i], 
         perc_digits+5, perc_digits, data[i]/data[0] );
   }
   fprintf(fp, "\n");

   // GEN gen_data[23];
   // gen_data[0]    = GEN_Create( &t_prof->name,         DATATYPE_STRING,  sizeof(char*) );
   // gen_data[1]    = GEN_Create( &q_seq->name,          DATATYPE_STRING,  sizeof(char*) );
   // gen_data[2]    = GEN_Create( &worker->result->time, DATATYPE_FLOAT,   sizeof(float) );
   // gen_data[3]    = GEN_Create( &times->load_target,   DATATYPE_FLOAT,   sizeof(float) );  
   // gen_data[4]    = GEN_Create( &times->load_query,    DATATYPE_FLOAT,   sizeof(float) );
   // gen_data[5]    = GEN_Create( &times->lin_cloud_fwd, DATATYPE_FLOAT,   sizeof(float) );
   // gen_data[6]    = GEN_Create( &times->lin_cloud_bck, DATATYPE_FLOAT,   sizeof(float) );
   // gen_data[7]    = GEN_Create( &times->lin_merge,     DATATYPE_FLOAT,   sizeof(float) );
   // gen_data[8]    = GEN_Create( &times->lin_reorient,  DATATYPE_FLOAT,   sizeof(float) );
   // gen_data[9]    = GEN_Create( &times->lin_bound_fwd, DATATYPE_FLOAT,   sizeof(float) );
   // gen_data[10]   = GEN_Create( &times->lin_bound_bck, DATATYPE_FLOAT,   sizeof(float) );
   // gen_data[11]   = GEN_Create( &times->sp_build_mx,   DATATYPE_FLOAT,   sizeof(float) );
   // gen_data[12]   = GEN_Create( &times->sp_bound_fwd,  DATATYPE_FLOAT,   sizeof(float) );
   // gen_data[13]   = GEN_Create( &times->sp_bound_bck,  DATATYPE_FLOAT,   sizeof(float) );
   // gen_data[14]   = GEN_Create( &times->sp_posterior,  DATATYPE_FLOAT,   sizeof(float) );
   // gen_data[15]   = GEN_Create( &times->sp_decodedom,  DATATYPE_FLOAT,   sizeof(float) );
   // gen_data[16]   = GEN_Create( &times->sp_biascorr,   DATATYPE_FLOAT,   sizeof(float) );
   // gen_data[17]   = GEN_Create( &times->sp_optacc,     DATATYPE_FLOAT,   sizeof(float) );
   // gen_data[18]   = GEN_Create( &times->dom_bound_fwd, DATATYPE_FLOAT,   sizeof(float) );
   // gen_data[19]   = GEN_Create( &times->dom_bound_bck, DATATYPE_FLOAT,   sizeof(float) );
   // gen_data[20]   = GEN_Create( &times->dom_posterior, DATATYPE_FLOAT,   sizeof(float) );
   // gen_data[21]   = GEN_Create( &times->dom_biascorr,  DATATYPE_FLOAT,   sizeof(float) );
   // gen_data[22]   = GEN_Create( &times->dom_optacc,    DATATYPE_FLOAT,   sizeof(float) );

   // REPORT_entry_multiline( fp, headers, gen_data, num_fields, sig_digits );
}