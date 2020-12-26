/*******************************************************************************
 *  FILE:      pipeline_mmseqs_plus.c
 *  PURPOSE:   Cloud Search Pipeline for MMSEQS: R.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
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
#include "../utilities/utilities.h"
#include "../objects/objects.h"
#include "../parsers/parsers.h"
#include "../algs_linear/algs_linear.h"
#include "../algs_quad/algs_quad.h"
#include "../algs_naive/algs_naive.h"

/* header */
#include "pipelines.h"

/* mmseqs pipeline */
void mmseqs_plus_pipeline( WORKER* worker )
{
   printf("# begin mmseqs pipeline...\n");
   
   /* command and options */
   char* command[27];
   for (int i = 0; i < 27; i++) {
      command[i] = ERROR_malloc( sizeof(char) * MAX_PATH_LEN );
      command[i][0] = '\0';
   }
   command[26] = NULL;
   /* NOTE: No filepath allowed exceed MAX_PATH_LEN */ 

   /* Variable Order */
   /* -- SCRIPTS --
    * (0)  MMORE-SEQS Script Location (executable)
    * -- MAIN ARGS --
    * (1)  Target File 
    * (2)  Query File
    * (3)  Target File Type 
    * (4)  Query File Type
    * -- OPTIONS --
    * (5)  Results Output File 
    * (6)  Remove Temp? 
    * (7)  Verbose Level
    * -- MMseqs Args -- 
    * (8)  Kmer Length
    * (9)  Kmer Score (Prefilter Threshold)
    * (10)  Minimum Ungapped Score (Ungapped Viterbi Threshold)
    * (11)  P-Value Cutoff (Gapped Viteri / MMseqs Reporting Threshold)
    * (12) E-Value Cutoff (Gapped Viterbi / MMseqs Reporting Threshold)
    * -- MMORE-SEQS / FB-PRUNER Args --
    * (13)  Alpha 
    * (14)  Beta
    * (15)  Gamma
    * (16)  P-Value Reporting Threshold
    * (17)  E-Value Reporting Threshold
    * (18)  Composition Bias
    * -- TOOL LOCATIONS --
    * (19)  FASTA-TO-HMM Script Location (executable)
    * (20)  MMSEQS Binary Location (executable)
    * (21)  HMMBUILD Binary Location (executable)
    * -- OUTPUTS --
    * (22)  Std-Output
    * (23)  Tbl-Output
    * (23)  M8-Output
    * (24)  MY-Output
    * -- NULL TERMINATED --
    * (25)  NULL
    */ 

   /* -- MAIN COMMAND LOCATION - */
   snprintf(command[0], MAX_PATH_LEN, " %s ", MMSEQS_PLUS_EASY_SCRIPT);

   /* --- MAIN ARGS --- */
   /* Target File Path */
   snprintf(command[1], MAX_PATH_LEN, " %s ", worker->args->q_filepath);
   /* Query File Path */
   snprintf(command[2], MAX_PATH_LEN, " %s ", worker->args->t_filepath);
   /* Target File Type */
   snprintf(command[3], MAX_PATH_LEN, " %d ", worker->args->q_filetype);
   /* Query File Type */
   snprintf(command[4], MAX_PATH_LEN, " %d ", worker->args->t_filetype);

   /* --- OPTIONS --- */
   /* Results File */
   snprintf(command[5], MAX_PATH_LEN, " %s ", worker->args->output_filepath);
   /* Remove Temporary Files */
   snprintf(command[6], MAX_PATH_LEN, " %d ", worker->args->tmp_remove);
   /* Verbosity Level */
   snprintf(command[7], MAX_PATH_LEN, " %d ", worker->args->verbose_level);

   /* --- MMSEQS ARGS --- */
   /* Kmer Length (Pre-filter) */
   snprintf(command[8], MAX_PATH_LEN, " %d ", worker->args->mmseqs_kmer);
   /* K-Score (Pre-filter) */
   snprintf(command[9], MAX_PATH_LEN, " %d ", worker->args->mmseqs_prefilter);
   /* Ungapped Viterbi */
   snprintf(command[10], MAX_PATH_LEN, " %d ", worker->args->mmseqs_ungapped_vit);
   /* E-Value MMseqs Reporting Threshold (Gapped Viterbi) */
   snprintf(command[11], MAX_PATH_LEN, " %f ", worker->args->mmseqs_evalue);
   /* P-Value MMseqs Reporting Threshold (Gapped Viterbi) */
   snprintf(command[12], MAX_PATH_LEN, " %f ", worker->args->mmseqs_pvalue);
  

   /* --- MMORE-SEQS / FB-PRUNER ARGS --- */
   /* Cloud Search alpha */
   snprintf(command[13], MAX_PATH_LEN, " %f ", worker->args->alpha);
   /* Cloud Search beta */
   snprintf(command[14], MAX_PATH_LEN, " %f ", worker->args->beta);
   /* Cloud Search gamma */
   snprintf(command[15], MAX_PATH_LEN, " %d ", worker->args->gamma);
   /* E-Value FB-Pruner Reporting Threshold */
   snprintf(command[16], MAX_PATH_LEN, " %f ", worker->args->mmseqs_evalue);
   /* P-Value FB-Pruner Reporting Threshold */
   snprintf(command[17], MAX_PATH_LEN, " %f ", worker->args->mmseqs_pvalue);  
   /* P-Value FB-Pruner Reporting Threshold */
   snprintf(command[18], MAX_PATH_LEN, " %d ", worker->args->is_compo_bias); 

   /* --- OUTPUT PATHS --- */
   /* standard output */
   snprintf(command[19], MAX_PATH_LEN, " %s ", worker->args->output_filepath);
   /* standard output */
   snprintf(command[20], MAX_PATH_LEN, " %s ", (worker->args->is_tblout ? worker->args->tblout_filepath : "--") );
   /* m8 output */
   snprintf(command[21], MAX_PATH_LEN, " %s ", (worker->args->is_m8out ? worker->args->m8out_filepath : "--") );
   /* my output */
   snprintf(command[22], MAX_PATH_LEN, " %s ", (worker->args->is_myout ? worker->args->myout_filepath : "--") );

   /* --- TOOL LOCATIONS --- */
   /* fasta-to-hmm script */
   snprintf(command[23], MAX_PATH_LEN, " %s ", FASTA_TO_HMM_SCRIPT);
   /* mmseqs binary */
   snprintf(command[24], MAX_PATH_LEN, " %s ", MMSEQS_BIN);
   /* hmmbuild binary */
   snprintf(command[25], MAX_PATH_LEN, " %s ", HMMBUILD_BIN);

   /* -- NULL TERMINATOR --- */
   command[26] = NULL;

   for (int i = 0; i < 27; i++) {
      printf("[%d] %s\n", i, command[i] );
   }
   
   /* execute script */
   // int exit_code = execvp( command[0], command );
   
   // char* test_cmd[] = { command[0], command[1], command[2], NULL };
   char* test_cmd[] = { MMSEQS_PLUS_EASY_SCRIPT, command[0] };
   int exit_code = execvp( command[0], command );

   /* code should not reach here */
   printf("# EXIT CODE = %d.  COMMAND FAILED.\n", exit_code);
   exit(EXIT_FAILURE);
}

