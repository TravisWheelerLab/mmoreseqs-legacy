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

/* local imports */
#include "../objects/structs.h"
#include "../utilities/_utilities.h"
#include "../objects/_objects.h"

/* header */
#include "_reporting.h"
#include "domtblout.h"

/* === DOMTBLOUT OUTPUT === */
/* EXAMPLE:
 *
      [header/]
      #                                                                   --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
      # target name            accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
      #    ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
      [entry/]
      3-PAP/16/510-647/718-827 -          3-PAP                PF12578.1    3.2e-12   34.5   0.0   1.5e-08   22.6   0.0   2.5   2   0   0   2   2   2   2 domains: MTMRA_DANRE/548-685 C3Z9W9_BRAFL/506-615
      3-PAP/13/1-136/374-513   -          3-PAP                PF12578.1    2.7e-11   31.5   0.0   3.9e-10   27.7   0.0   2.4   2   0   0   2   2   2   1 domains: MTMRB_MOUSE/553-688 B7Q8P0_IXOSC/504-643
      3-PAP/14/86-218/365-501  -          3-PAP                PF12578.1    1.1e-07   19.9   0.0   2.1e-07   18.9   0.0   1.5   1   0   0   1   1   1   1 domains: MTMRC_PONAB/559-691 A4HUS9_LEIIN/8-144
      [footer/]
      #
      # Program:         hmmsearch
      # Version:         3.3 (Nov 2019)
      # Pipeline mode:   SEARCH
      # Query file:      test-input/3-PAP.hmm
      # Target file:     test-input/3-PAP.fa
      # Option settings: hmmsearch --tblout tblout.csv test-input/3-PAP.hmm test-input/3-PAP.fa 
      # Current dir:     /home/devreckas/Google-Drive/Wheeler-Labs/Personal_Work/fb-pruner/fb-pruner
      # Date:            Fri Aug 28 00:42:33 2020
      # [ok]
 *
 */

/*!   FUNCTION:   REPORT_domtblout_header()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
STATUS_FLAG 
REPORT_domtblout_header(   WORKER*  worker,
                           FILE*    fp )
{
   int qnamew = 20;
   int tnamew = 20;
   int qaccw  = 10;
   int taccw  = 10;
   int posw   = 7;

   fprintf( fp, "#%*s %22s %22s %33s\n",
      tnamew + qnamew + taccw + qaccw + 2, 
      "", 
      "--- full sequence ----", 
      "--- best 1 domain ----", 
      "--- domain number estimation ----" 
   );

   fprintf( fp, "#%-*s %-*s %-*s %-*s %9s %6s %5s %9s %6s %5s %5s %3s %3s %3s %3s %3s %3s %3s %s\n",
      tnamew - 1, 
      " target name",        
      taccw, 
      "accession",  
      qnamew, 
      "query name",           
      qaccw, 
      "accession",  
      "  E-value", 
      " score", 
      " bias", 
      "  E-value", 
      " score", 
      " bias", 
      "exp", 
      "reg", 
      "clu", 
      " ov", 
      "env", 
      "dom", 
      "rep", 
      "inc", 
      "description of target" 
   );

   fprintf( fp, "#%*s %*s %*s %*s %9s %6s %5s %9s %6s %5s %5s %3s %3s %3s %3s %3s %3s %3s %s\n",
      tnamew - 1, 
      "-------------------", 
      taccw, 
      "----------", 
      qnamew, 
      "--------------------", 
      qaccw, 
      "----------", 
      "---------", 
      "------", 
      "-----", 
      "---------", 
      "------", 
      "-----", 
      "---", 
      "---", 
      "---", 
      "---", 
      "---", 
      "---", 
      "---", 
      "---", 
      "---------------------" 
   );
}

/*!   FUNCTION:   REPORT_domtblout_entry()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
STATUS_FLAG 
REPORT_domtblout_entry(    WORKER*  worker,
                           RESULT*  result,
                           FILE*    fp )
{
   int qnamew = 20;
   int tnamew = 20;
   int qaccw  = 10;
   int taccw  = 10;
   int posw   = 7;

   /* computed output */
   float bias_correction = result->final_scores.pre_sc - result->final_scores.seq_sc;

   fprintf( fp, "%-*s %-*s %-*s %-*s %9.2g %6.1f %5.1f %9.2g %6.1f %5.1f %5.1f %3d %3d %3d %3d %3d %3d %3d %s\n", 
      /* query / target data */
      qnamew, worker->q_seq->name,                          /* query name */
      qaccw,  (NULL ? worker->q_seq->acc : "-"),            /* query accession */
      tnamew, worker->t_prof->name,                         /* target name */
      taccw,  (NULL ? worker->t_prof->acc : "-"),           /* target accession */
      /* full sequence */
      result->final_scores.eval,                            /* evalue */
      result->final_scores.nat_sc,                          /* score */
      bias_correction,                                      /* bias correction (pre_score - score) */
      /* best domain */
      0.0,                                                  /* evalue */
      0.0,                                                  /* score */
      0.0,                                                  /* bias correction */
      /* domain number estimation */
      0.0,                                                  /* number expected */
      0,                                                    /* number regions */
      0,                                                    /* number clustered */
      0,                                                    /* number overlapped */
      0,                                                    /* number envelopes */
      0,                                                    /* number domains */
      0,                                                    /* number reported */
      0,                                                    /* number included */
      /* target description */
      (NULL ? worker->t_prof->desc : "--")                  /* query description */
   );
}

/*!   FUNCTION:   REPORT_domtblout_footer()
 *    SYNOPSIS:   Print footer
 *                (modeled after HMMER, see example)              
 */
STATUS_FLAG 
REPORT_domtblout_footer(   WORKER*  worker,
                           FILE*    fp )
{
   ARGS* args     = worker->args;
   int   left_pad = 20;

   /* get current directory */
   char cwd[256];
   if ( getcwd(cwd, sizeof(cwd)) == NULL ) return STATUS_FAILURE;
   /* get current date/time */
   char* time_str = NULL;
   // time_str = CLOCK_GetDateTimeString( NULL );

   fprintf( fp, "# \n");
   fprintf( fp, "# %*s %s\n",       left_pad,   "Program:",          BUILD_PROGRAM);
   fprintf( fp, "# %*s %s (%s)\n",  left_pad,   "Version:",          BUILD_VERSION, BUILD_DATE );
   fprintf( fp, "# %*s %s\n",       left_pad,   "Pipeline mode:",    PIPELINES[args->pipeline_mode].name );
   fprintf( fp, "# %*s %s\n",       left_pad,   "Query file:",       args->q_filein );
   fprintf( fp, "# %*s %s\n",       left_pad,   "Target file:",      args->t_filein );
   fprintf( fp, "# %*s %s\n",       left_pad,   "Option settings:",  "" );
   if (cwd) {
   fprintf( fp, "# %*s %s\n",       left_pad,   "Current dir:",      cwd );
   }
   // if (time_str)  fprintf( fp, "# %*s %s\n",       left_pad,   "Date:",             time_str );
   fprintf( fp, "# %*s\n",          left_pad,   "[ok]" );

   return STATUS_SUCCESS;
}

