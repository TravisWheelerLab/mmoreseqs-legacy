/*******************************************************************************
 *     FILE:   arg_parser.h
 *  PURPOSE:   Parses command line arguments. 
 *
 *   AUTHOR:   Dave Rich
 *     BUGS:
 *       - None Known.
 *    NOTES:
 *       - Working on making argument options a generic datatype rather than an if-else ladder.
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

/* header */
#include "_parsers.h"
#include "arg_parser.h"

/* private functions */

void
ARGS_MainArg_Parser(    ARGS*          args,
                        COMMANDLINE*   cmd );


/*! FUNCTION:  ARGS_Parse()
 *  SYNOPSIS:  Parses Arguments from the command line.
 */
void   
ARGS_Parse(    ARGS*          args,
               int            argc, 
               char*          argv[],
               COMMANDLINE*   cmd,
               ARG_OPTS*      arg_opts )
{
   int   num_main_args     = 2; 
   char* flag              = NULL;
   /* required arguments */
   int   req_args          = 0;
   /* remaining arguments counter */
   int   args_rem          = argc - 1;
   /* current argument index */
   int   arg_cur           = 1;

   printf("NUM_ARGS: %d\n", argc);

   ARGS_SetDefaults( args );
   ARGS_SetOptions( args, arg_opts );

   /* if no <command> argument given, run test case if in debug mode */
   if (argc <= 1) {
      printf("Usage: mmore <command> <main_args...>\n");
      printf("Hint: For more information, use '-h'");
      ERRORCHECK_exit(EXIT_FAILURE);
   }

   /* check for help flag */
   if ( STR_Compare( argv[1], "-h") == 0 || STR_Compare( argv[1], "--help") == 0 ) {
      ARGS_Help_Info();
   }

   /* first argument is command pipeline (id_number or name) */
   bool found_pipeline = false;
   if ( isdigit(argv[1][0]) ) {
      args->pipeline_mode = atoi(argv[1]);
   } else {
      for (int i = 0; i < NUM_PIPELINE_MODES; i++) {
         PIPELINE* pipeline = &PIPELINES[i];
         if ( STR_Compare( argv[1], pipeline->name ) == 0 ) {
            args->pipeline_mode = i;
            args->pipeline_name = STR_Create( pipeline->name );
            found_pipeline = true;
            break;
         }
      }
   }
   /* check that valid pipeline mode was entered */
   if ( found_pipeline == false || (args->pipeline_mode < 0) || (args->pipeline_mode >= NUM_PIPELINE_MODES) ) {
      fprintf(stderr, "ERROR: invalid pipeline/command was given: (%s, %d).\n", argv[1], args->pipeline_mode);
      fprintf(stderr, "VALID PIPELINE/COMMANDS OPTS: [ ");
      for (int i = 0; i < NUM_PIPELINE_MODES; i++) 
         fprintf(stderr, "%s, ", PIPELINES[i].name);
      fprintf(stderr, "]\n");
      ERRORCHECK_exit(EXIT_FAILURE);
   }
   args_rem    -= 1;
   arg_cur     += 1;

   /* set number of main arguments based on given pipeline */
   num_main_args = PIPELINES[args->pipeline_mode].num_main_args;
   /* check proper number of main args remain */
   if ( args_rem < num_main_args ) {
      fprintf(stderr, "ERROR: Improper number of main args. [required: %d/%d]\n", args_rem, num_main_args);
      #if DEBUG 
         printf("Using DEFAULT arguments...\n\n");
         return;
      #else 
         ERRORCHECK_exit(EXIT_FAILURE);
      #endif
   }

   /* check if proper number of args */
   for (int i = 0; i < num_main_args; i++) {
      if ( STR_StartsWith(argv[i + 1], "--") == 0 ) {
         fprintf(stderr, "ERROR: Improper number of main arguments.\n");
         ERRORCHECK_exit(EXIT_FAILURE);
      }
   }

   /* parse main commands */
   if   ( STR_Equals( args->pipeline_name, "search" ) )
   {
      args->t_mmore_filein     = STR_Set( args->t_mmore_filein,     argv[2] );
      args->q_mmore_filein     = STR_Set( args->q_mmore_filein,     argv[3] );
      args->t_mmseqs_p_filein  = STR_Set( args->t_mmseqs_p_filein,  argv[4] );
      args->t_mmseqs_s_filein  = STR_Set( args->t_mmseqs_s_filein,  argv[5] );
      args->q_mmseqs_filein    = STR_Set( args->q_mmseqs_filein,    argv[6] );

      args->t_mmore_p_filetype   = FILE_HMM;
      args->q_filetype           = FILE_FASTA;
      args->t_mmseqs_p_filetype  = FILE_MMDB_P;
      args->t_mmseqs_s_filetype  = FILE_MMDB_S;
      args->q_mmseqs_filetype    = FILE_MMDB_S;
   }
   elif ( STR_Equals( args->pipeline_name, "mmore-search" ) )
   {
      args->t_filein           = STR_Set( args->t_filein,           argv[2] );
      args->q_filein           = STR_Set( args->q_filein,           argv[3] );
      args->mmseqs_m8_filein   = STR_Set( args->mmseqs_m8_filein,   argv[4] );

      args->t_filetype = FILE_HMM;
      args->q_filetype = FILE_FASTA;
   }
   elif ( STR_Equals( args->pipeline_name, "mmseqs-search" ) )
   {
      args->t_mmseqs_p_filein  = STR_Set( args->t_mmseqs_p_filein,  argv[2] );
      args->t_mmseqs_s_filein  = STR_Set( args->t_mmseqs_s_filein,  argv[3] );
      args->q_mmseqs_filein    = STR_Set( args->q_mmseqs_filein,    argv[4] );

      args->t_mmseqs_p_filetype  = FILE_MMDB_P;
      args->t_mmseqs_s_filetype  = FILE_MMDB_S;
      args->q_mmseqs_filetype    = FILE_MMDB_S;
   }
   elif ( STR_Equals( args->pipeline_name, "prep" ) )
   {
      args->target_prep       = STR_Set( args->target_prep,       argv[2] );
      args->query_prep        = STR_Set( args->query_prep,        argv[3] );
      args->tmp_folderpath    = STR_Set( args->tmp_folderpath,    argv[4] );
      args->prep_folderpath   = STR_Set( args->prep_folderpath,   argv[4] );

      args->target_prep_type  = FILE_MSA;
      args->query_prep_type   = FILE_FASTA;
   }
   elif ( STR_Equals( args->pipeline_name, "prep-search" ) )
   {
      
      args->prep_folderpath   = STR_Set( args->prep_folderpath,   argv[2] );
      args->tmp_folderpath    = STR_Set( args->tmp_folderpath,    argv[2] );
   } 
   elif ( STR_Equals( args->pipeline_name, "easy-search" ) )
   {
      args->target_prep       = STR_Set( args->target_prep,       argv[2] );
      args->query_prep        = STR_Set( args->query_prep,        argv[3] );
      args->tmp_folderpath    = STR_Set( args->tmp_folderpath,    argv[4] );
      args->prep_folderpath   = STR_Set( args->prep_folderpath,   argv[4] );

      args->target_prep_type  = FILE_MSA;
      args->query_prep_type   = FILE_FASTA;
   }
   elif ( STR_Equals( args->pipeline_name, "index" ) )
   {
      args->t_index_filein    = STR_Set( args->t_index_filein,    argv[2] );
      // args->t_indexout     = STR_Set( args->t_indexout,     argv[3] );
   }
   else {
      fprintf(stderr, "ERROR: Pipeline option '%s' is currently not supported.\n", args->pipeline_name);
      ERRORCHECK_exit(EXIT_FAILURE);
   }
   args_rem    -= num_main_args;
   arg_cur     += num_main_args;

   /* parse remaining flags and options */
   for (int i = arg_cur; i < argc; ++i)
   {
      /* if long flag */
      if ( STR_ComparePrefix( argv[i], "--", 2 ) == 0 ) 
      {
         /* === PIPELINE OPTIONS === */
         if      ( STR_Compare( argv[i], (flag = "--help") ) == 0 ) {
            ARGS_Help_Info();
         }
         else if ( STR_Compare( argv[i], (flag = "--verbose") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->verbose_level = atoi(argv[i]);
               if ( args->verbose_level < 0 || args->verbose_level > NUM_VERBOSITY_MODES - 1 ) {
                  fprintf(stderr, "ERROR: Verbose Level (%d) is outside acceptable range (%d,%d).\n", 
                     args->verbose_level, 0, NUM_VERBOSITY_MODES - 1 );
                  args->verbose_level = MAX( args->verbose_level, 0 );
                  args->verbose_level = MIN( args->verbose_level, NUM_VERBOSITY_MODES-1 );
                  fprintf(stderr, "WARNING: Verbose Level set to: %d.\n", args->verbose_level );
               } 
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--num-threads") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->num_threads = atoi(argv[i]);
               if ( args->num_threads < 1 ) {
                  fprintf(stderr, "ERROR: Number of Threads (%d) is outside acceptable range (%d,%d).\n", 
                     args->num_threads, 1, 16 );
                  args->num_threads = MAX( args->num_threads, 1 );
                  args->num_threads = MIN( args->num_threads, 16 );
                  fprintf(stderr, "WARNING: Number of Threads set to: %d.\n", args->num_threads );
               } 
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--enforce-errors") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               if ( atoi(argv[i]) == 0 ) {
                  args->enforce_warnings = false;
               } else if ( atoi(argv[i]) == 1 ) {
                  args->enforce_warnings = true;
               }
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--enforce-warnings") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               if ( atoi(argv[i]) == 0 ) {
                  args->enforce_warnings = false;
               } else if ( atoi(argv[i]) == 1 ) {
                  args->enforce_warnings = true;
               }
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--search-type") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               if   ( STR_Equals( argv[i], "P2S" ) == true ) {
                  args->search_name       = STR_Set( args->search_name, argv[i] );
                  args->target_prep_type  = FILE_MSA; 
                  args->query_prep_type   = FILE_FASTA;
               }
               elif ( STR_Equals( argv[i], "S2S" ) == true ) {
                  args->search_name       = STR_Set( args->search_name, argv[i] );
                  args->target_prep_type  = FILE_FASTA;
                  args->query_prep_type   = FILE_FASTA;
               } 
               else {
                  fprintf(stderr, "ERROR: '%s' is not a valid argument for --search-type.\n", argv[i] );
                  ERRORCHECK_exit(EXIT_FAILURE);
               }
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         /* === INPUT FILES === */
         else if ( STR_Compare( argv[i], (flag = "--index") ) == 0 ) {
            req_args = 2;
            if (i+req_args < argc) {
               i++;
               args->t_index_filein = STR_Set( args->t_index_filein, argv[i] );
               i++;
               args->q_index_filein = STR_Set( args->q_index_filein, argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument=.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--local-tools") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->is_use_local_tools = atoi(argv[i]); 
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument=.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--guess-ftype") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->is_guess_filetype = atoi(argv[i]); 
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument=.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--mmoreseqs-ftype") ) == 0 ) {
            req_args = 5;
            if (i+req_args < argc) {
               i++;
               args->t_filetype = atoi(argv[i]); 
               i++;
               args->q_filetype = atoi(argv[i]);
               i++;
               args->t_mmseqs_p_filetype = atoi(argv[i]);
                i++;
               args->t_mmseqs_s_filetype = atoi(argv[i]);
               i++;
               args->q_mmseqs_filetype = atoi(argv[i]);
               /* since filetype supplied, turn off guesser */
               args->is_guess_filetype = false;
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--mmoreseqs-main-ftype") ) == 0 ) {
            req_args = 2;
            if (i+req_args < argc) {
               i++;
               args->t_filetype = atoi(argv[i]); 
               i++;
               args->q_filetype = atoi(argv[i]);
               /* since filetype supplied, turn off guesser */
               args->is_guess_filetype = false;
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--tmp") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->tmp_folderpath = STR_Set( args->tmp_folderpath, argv[i] );
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--prep") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->tmp_folderpath = STR_Set( args->prep_folderpath, argv[i] );
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--mmseqs-m8") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->mmseqs_m8_filein = STR_Set( args->mmseqs_m8_filein, argv[i] );
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         /* === INPUT DATA === */
         else if ( STR_Compare( argv[i], (flag = "--dbsizes") ) == 0 ) {
            req_args = 2;
            if (i+req_args < argc) {
               i++;
               args->q_dbsize = atoi(argv[i]);
               i++;
               args->t_dbsize = atoi(argv[i]);
            } else {
               fprintf(stderr, "ERROR: '%s' flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         /* === MMORE PARAMETERS === */
         else if ( STR_Compare( argv[i], (flag = "--alpha") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->alpha = atof(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--beta") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->beta = atof(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--gamma") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->gamma = atoi(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--mmoreseqs-viterbi-pval") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->threshold_vit = atof(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--mmoreseqs-cloud-pval") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->threshold_cloud = atof(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--mmoreseqs-boundfwd-pval") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->threshold_boundfwd = atof(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--mmoreseqs-fwd-pval") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->threshold_fwd = atoi(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--eval") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->mmseqs_evalue = atof(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         /* ==== MMORE OPTIONS === */
         else if ( STR_Compare( argv[i], (flag = "--run-prep") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->tmp_folderpath = STR_Set( args->tmp_folderpath, argv[i] );
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--run-filter") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->is_run_filter = atoi(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--run-bias") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->is_run_bias = atoi(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--run-full") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->is_run_full = atoi(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--run-domains") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->is_run_domains = atoi(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--run-mmseqsaln") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->is_run_mmseqsaln = atoi(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--run-vitaln") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->is_run_vitaln = atoi(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--run-postaln") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->is_run_postaln = atoi(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--run-mmseqs-ungapped") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->is_run_mmseqs_ungapped = atoi(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--run-mmore") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->is_run_mmore = atoi(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         /* === MMSEQS PARAMETERS === */
         else if ( STR_Compare( argv[i], (flag = "--mmseqs-split") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->mmseqs_split = atoi(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--mmseqs-kmer") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->mmseqs_kmer = atoi(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--mmseqs-ungapped-vit") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->mmseqs_ungapped_vit = atoi(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--mmseqs-eval") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->mmseqs_evalue = atof(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--mmseqs-pval") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->mmseqs_pvalue = atof(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--mmseqs-hits-per-search") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->mmseqs_maxhits = atoi(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         /* === MMSEQS DATA === */
         else if ( STR_Compare( argv[i], (flag = "--mmseqs-times") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->mmseqs_maxhits = atoi(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         /* === SEARCH/RANGE OPTIONS === */
         else if ( STR_Compare( argv[i], (flag = "--range") ) == 0 ) {
            req_args = 2;
            if (i+req_args < argc) {
               i++;
               args->list_range.beg = atoi(argv[i]);
               i++;
               args->list_range.end = atoi(argv[i]);
            } else {
               fprintf(stderr, "ERROR: '%s' flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--search-mode") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->search_mode = atoi(argv[i]);
               if ( args->search_mode < 0 || args->search_mode > NUM_SELECT_SEARCHES ) {
                  fprintf(stderr, "ERROR: Verbose level (%d) is outside acceptable range (%d,%d).\n", 
                     args->search_mode, 0, NUM_SELECT_SEARCHES );
                  args->search_mode = MAX( args->search_mode, 0 );
                  args->search_mode = MIN( args->search_mode, NUM_SELECT_SEARCHES-1 );
               } 
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         /* === INTERRIM OUTPUT === */
         else if ( STR_Compare( argv[i], (flag = "--mmseqs-p2sout") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               free(args->stdout_fileout);
               args->stdout_fileout = STR_Create(argv[i]);
               args->is_redirect_stdout = true;
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         /* === OUTPUT === */
         else if ( STR_Compare( argv[i], (flag = "--stderr") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               fclose(stdout);
               free(args->stderr_fileout);
               args->stderr_fileout = STR_Create(argv[i]);
               args->is_redirect_stderr = true;
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--stdout") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               free(args->stdout_fileout);
               args->stdout_fileout = STR_Create(argv[i]);
               args->is_redirect_stdout = true;
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--allout") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               STR_Destroy( args->hmmerout_fileout );
               args->hmmerout_fileout = STR_Concat(argv[i], ".hmmerout");
               STR_Destroy( args->m8out_fileout );
               args->m8out_fileout = STR_Concat(argv[i], ".m8out");
               args->is_m8out = true;
               STR_Destroy(args->myout_fileout);
               args->myout_fileout = STR_Concat(argv[i], ".myout");
               args->is_myout = true;
               STR_Destroy(args->mydom_fileout);
               args->mydom_fileout = STR_Concat(argv[i], ".mydomout");
               args->is_mydom = true;
               STR_Destroy(args->mytime_fileout);
               args->mytime_fileout = STR_Concat(argv[i], ".mytimeout");
               args->is_mytimeout = true;
               STR_Destroy(args->mythresh_fileout);
               args->mythresh_fileout = STR_Concat(argv[i], ".mythreshout");
               args->is_mythreshout = true;
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--domtblout") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               ERROR_free(args->hmmerout_fileout);
               args->hmmerout_fileout = STR_Create(argv[i]);
               args->is_hmmerout = true;
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--m8out") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               ERROR_free(args->m8out_fileout);
               args->m8out_fileout = STR_Create(argv[i]);
               args->is_m8out = true;
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--myout") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               ERROR_free(args->myout_fileout);
               args->myout_fileout = STR_Create(argv[i]);
               args->is_myout = true;
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--mydomtblout") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               ERROR_free(args->mydom_fileout);
               args->mydom_fileout = STR_Create(argv[i]);
               args->is_mydom = true;
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--mytimeout") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               ERROR_free(args->mytime_fileout);
               args->mytime_fileout = STR_Create(argv[i]);
               args->is_mytimeout = true;
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--mythreshout") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               ERROR_free(args->mythresh_fileout);
               args->mythresh_fileout = STR_Create(argv[i]);
               args->is_mythreshout = true;
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--output") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               args->stdout_fileout = STR_Create(argv[i]);
               args->is_redirect_stdout = true;
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--customout") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               ERROR_free(args->hmmerout_fileout);
               args->hmmerout_fileout = STR_Create(argv[i]);
               args->is_hmmerout = true;
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--debugout") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               ERROR_free(args->dbg_folderpath);
               args->dbg_folderpath = STR_Create(argv[i]);
               args->is_debug = true;
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         /* === NOT PROPER FLAG === */
         else {
            fprintf(stderr, "ERROR: '%s' is not a recognized flag.\n", argv[i]);
            ERRORCHECK_exit(EXIT_FAILURE);
         }
      }
      else
      {
         fprintf(stderr, "ERROR: '%s' is not associated with an argument or flag\n", argv[i]);
         ERRORCHECK_exit(EXIT_FAILURE);
      }
   }

   /* print commandline input */
   
}

/*! FUNCTION:  ARGS_ArgOpt_Parser()
 *  SYNOPSIS:  Set default arguments.
 */
void
ARGS_ArgOpt_Parser(  ARGS*          args, 
                     COMMANDLINE*   cmd, 
                     ARG_OPTS*      argopt )
{

}

/*! FUNCTION:  ARGS_MainArg_Parser()
 *  SYNOPSIS:  Set default arguments.
 */
void
ARGS_MainArg_Parser(    ARGS*          args,
                        COMMANDLINE*   cmd )
{

}

/*! FUNCTION:  ARGS_SetDefaults()
 *  SYNOPSIS:  Set default arguments.
 */
void  
ARGS_SetDefaults( ARGS*    args )
{
   /* --- COMMANDLINE --- */

   /* --- PIPELINE MODE --- */
   args->pipeline_mode           = PIPELINE_TEST;
   args->pipeline_name           = NULL;
   args->verbose_level           = VERBOSE_LOW;
   args->num_threads             = 1;
   args->search_mode             = MODE_UNILOCAL;
   args->qt_search_space         = SELECT_ALL_V_ALL;
   args->tmp_folderpath          = NULL;
   args->tmp_remove              = false;

   /* --- PREPARATION OPTIONS --- */
   /* files/folders */
   args->is_run_prep             = true;
   args->is_prep_copy            = true;
   args->prep_folderpath         = NULL;
   args->target_prep             = NULL;
   args->query_prep              = NULL;
   /* types */
   args->target_prep_type        = FILE_MSA; 
   args->query_prep_type         = FILE_FASTA;

   /* --- PIPELINE OPTIONS --- */
   args->search_name             = STR_Create("P2S");
   args->is_run_bias             = true;
   args->is_run_pruned           = true;
   args->is_run_full             = false;
   args->is_run_domains          = true;
   args->is_run_mmseqsaln        = false;
   args->is_run_vitaln           = true;
   args->is_run_postaln          = false;
   args->is_run_mmseqs_ungapped  = false;
   args->is_run_mmore            = true;

   /* --- DEBUG OPTIONS --- */
   args->is_use_local_tools      = false;
   args->is_recycle_mx           = false;
   args->is_debug                = true;
   args->dbg_folderpath          = STR_Create(DEBUG_FOLDER "/");

   /* --- INPUT --- */
   /* filepath */
   args->t_filein              = NULL;
   args->q_filein              = NULL;
   args->t_mmore_filein        = NULL;
   args->q_mmore_filein        = NULL;
   args->t_mmseqs_p_filein     = NULL;
   args->t_mmseqs_s_filein     = NULL;
   args->q_mmseqs_filein       = NULL;
   /* filetype */
   args->is_guess_filetype       = true;
   args->t_filetype              = FILE_HMM;
   args->q_filetype              = FILE_FASTA;
   args->t_mmore_filetype        = FILE_HMM; 
   args->q_mmore_filetype        = FILE_FASTA;
   args->t_mmseqs_p_filetype     = FILE_MMDB_P;
   args->t_mmseqs_s_filetype     = FILE_MMDB_S;
   args->q_mmseqs_filetype       = FILE_MMDB_S;
   /* indexes */
   args->t_index_filein             = NULL;
   args->q_index_filein             = NULL;

   /* --- INTERRIM OUTPUT --- */
   args->mmseqs_m8_filein     = STR_Create("mmseqs.results.m8out");

   /* --- OUTPUT --- */
   /* default outputs */
   args->is_redirect_stdout      = false;
   args->stdout_fileout          = STR_Create("mmore.results.stdout");
   args->is_redirect_stderr      = false;
   args->stderr_fileout          = STR_Create("mmore.results.stderr");
   /* special outputs */
   args->is_allout               = false;
   args->allout_fileout          = STR_Create("mmore.results");
   args->is_hmmerout               = false;
   args->hmmerout_fileout          = STR_Create("mmore.results.tblout");
   args->is_m8out                = true;
   args->m8out_fileout           = STR_Create("mmore.results.m8");
   args->is_myout                = false;
   args->myout_fileout           = STR_Create("mmore.results.myout");
   args->is_mydom                = false;
   args->mydom_fileout           = STR_Create("mmore.results.mydomout");
   args->is_mytimeout            = false;
   args->mytime_fileout          = STR_Create("mmore.results.mytimeout");
   args->is_mythreshout          = false;
   args->mythresh_fileout        = STR_Create("mmore.results.mythreshout");
   // args->is_customout            = false;
   // args->customout_fileout      = STR_Create("results.customout");

   /* --- RANGE OPTIONS --- */
   args->t_range                 = (RANGE) { -1, -1 };    
   args->q_range                 = (RANGE) { -1, -1 };
   args->list_range              = (RANGE) { -1, -1 };

   /* --- MMORE / FB-PRUNER --- */
   args->alpha                   = 12.0f;
   args->beta                    = 20.0f;
   args->gamma                   = 5;
   args->mmore_evalue            = 2e2f;
   args->mmore_pvalue            = 1e-5f;

   /* --- MMSEQS --- */
   args->mmseqs_maxhits          = 1000;
   args->mmseqs_altalis          = 5;
   args->mmseqs_kmer             = 7;
   args->mmseqs_prefilter        = 80;
   args->mmseqs_ungapped_vit     = 0;
   args->mmseqs_evalue           = 1000.0;
   args->mmseqs_pvalue           = 1e-3f;

   /* --- SCORE FILTERS (P_VALUES) --- */
   args->is_run_filter           = false;
   args->threshold_pre           = 80;
   args->threshold_ungapped      = 0;
   args->threshold_p2s           = 1e-3f;
   args->threshold_s2s           = 1e-1f;
   args->threshold_vit           = 1e-3f;
   args->threshold_cloud         = 1e-5f;
   args->threshold_boundfwd      = 1e-5f;
   args->threshold_fwd           = 1e-5f;
}

/*! FUNCTION:  ARGS_Dump()
 *  SYNOPSIS:  Output arguments to <fp>.
 */
void 
ARGS_Dump(     ARGS*    args,
               FILE*    fp )
{
   int      pad               = 20;
   bool     align             = -1;     /* -1 for right alignment, 1 for left alignment */

   fprintf( fp, "# === MMORE OPTIONS =======================\n");
   /* --- PIPELINE --- */
   fprintf( fp, "# === PIPELINE ===\n");
   fprintf( fp, "# %*s:\t%s [%d]\n",    align * pad,  "PIPELINE",               PIPELINES[args->pipeline_mode].name,   args->pipeline_mode );
   fprintf( fp, "# %*s:\t%s [%d]\n",    align * pad,  "VERBOSITY_MODE",         VERBOSITY_NAMES[args->verbose_level],  args->verbose_level );
   fprintf( fp, "# === SCRIPTS ===\n");
   fprintf( fp, "# %*s:\t%s\n",         align * pad,  "PREP_SCRIPT",            PREP_SCRIPT );
   fprintf( fp, "# %*s:\t%s\n",         align * pad,  "PREPSEARCH_SCRIPT",      PREPSEARCH_SCRIPT );
   fprintf( fp, "# %*s:\t%s\n",         align * pad,  "SEARCH_SCRIPT",          SEARCH_SCRIPT );
   fprintf( fp, "# %*s:\t%s\n",         align * pad,  "EASYSEARCH_SCRIPT",      EASYSEARCH_SCRIPT );
   fprintf( fp, "# === PIPELINE OPTIONS ===\n");

   /* --- INPUT --- */
   fprintf( fp, "# === INPUT ===\n");
   fprintf( fp, "# %*s:\t%s [%s]\n",     align * pad,  "TARGET",                args->t_filein,             FILETYPE_NAME_Get( args->t_filetype ) );
   fprintf( fp, "# %*s:\t%s [%s]\n",     align * pad,  "QUERY",                 args->q_filein,             FILETYPE_NAME_Get( args->q_filetype ) );
   fprintf( fp, "# %*s:\t%s [%s]\n",     align * pad,  "TARGET_PREP",           args->target_prep,            FILETYPE_NAME_Get( args->target_prep_type ) );
   fprintf( fp, "# %*s:\t%s [%s]\n",     align * pad,  "QUERY_PREP",            args->query_prep,             FILETYPE_NAME_Get( args->query_prep_type ) );
   fprintf( fp, "# %*s:\t%s [%s]\n",     align * pad,  "TARGET_MMORE",          args->t_mmore_filein,       FILETYPE_NAME_Get( args->t_mmore_filetype ) );
   fprintf( fp, "# %*s:\t%s [%s]\n",     align * pad,  "QUERY_MMORE",           args->q_mmore_filein,       FILETYPE_NAME_Get( args->q_mmore_filetype ) );
   fprintf( fp, "# %*s:\t%s [%s]\n",     align * pad,  "TARGET_MMSEQS_P",       args->t_mmseqs_p_filein,    FILETYPE_NAME_Get( args->t_mmseqs_p_filetype ) );
   fprintf( fp, "# %*s:\t%s [%s]\n",     align * pad,  "TARGET_MMSEQS_S",       args->t_mmseqs_s_filein,    FILETYPE_NAME_Get( args->t_mmseqs_s_filetype ) );
   fprintf( fp, "# %*s:\t%s [%s]\n",     align * pad,  "QUERY_MMSEQS",          args->q_mmseqs_filein,      FILETYPE_NAME_Get( args->t_mmseqs_s_filetype ) );
   fprintf( fp, "# %*s:\t%s\n",          align * pad,  "T_INDEX_PATH",          args->t_index_filein );
   fprintf( fp, "# %*s:\t%s\n",          align * pad,  "Q_INDEX_PATH",          args->q_index_filein );
   fprintf( fp, "# %*s:\t%s\n",          align * pad,  "MMSEQS_M8",             args->mmseqs_m8_filein );
   fprintf( fp, "# %*s:\t%s\n",          align * pad,  "D2MP_FOLDER",            args->tmp_folderpath ); 
   fprintf( fp, "# %*s:\t%s\n",          align * pad,  "PREP_FOLDER",           args->prep_folderpath ); 
   fprintf( fp, "# \n" );
   /* --- MMSEQS --- */
   fprintf( fp, "# === MMSEQS ===\n");
   fprintf( fp, "# %*s:\t%d\n",          align * pad,  "MMSEQS_KMER",           args->mmseqs_kmer );
   fprintf( fp, "# %*s:\t%d\n",          align * pad,  "MMSEQS_KSCORE",         args->mmseqs_prefilter );
   fprintf( fp, "# %*s:\t%d\n",          align * pad,  "MMSEQS_UNGAPPED",       args->mmseqs_ungapped_vit );
   fprintf( fp, "# %*s:\t%.2e\n",        align * pad,  "MMSEQS_P2S_PVAL",       args->threshold_p2s );
   fprintf( fp, "# %*s:\t%.2e\n",        align * pad,  "MMSEQS_S2S_PVAL",       args->threshold_s2s ); 
   /* --- MMORE / FB-PRUNER --- */ 
   fprintf( fp, "# === MMORE ===\n");
   fprintf( fp, "# %*s:\t%.2f\n",        align * pad,  "MMORE_ALPHA",           args->alpha );
   fprintf( fp, "# %*s:\t%.2f\n",        align * pad,  "MMORE_BETA",            args->beta );
   fprintf( fp, "# %*s:\t%d\n",          align * pad,  "MMORE_GAMMA",           args->gamma );
   fprintf( fp, "# %*s:\t%.2e\n",        align * pad,  "MMORE_VITERBI_PVAL",    args->threshold_vit );
   fprintf( fp, "# %*s:\t%.2e\n",        align * pad,  "MMORE_CLOUD_PVAL",      args->threshold_cloud );
   fprintf( fp, "# %*s:\t%.2e\n",        align * pad,  "MMORE_BOUNDFWD_PVAL",   args->threshold_boundfwd );
   fprintf( fp, "# \n" );
   /* --- OUTPUT --- */
   fprintf( fp, "# === OUTPUT ===\n");
   fprintf( fp, "# %*s:\t%s\n",        align * pad,  "OUTPUT_FILEPATH",       args->stdout_fileout );
   fprintf( fp, "# %*s:\t%s\n",        align * pad,  "TBLOUT_FILEPATH",       args->hmmerout_fileout );
   fprintf( fp, "# %*s:\t%s\n",        align * pad,  "M8OUT_FILEPATH",        args->m8out_fileout );
   fprintf( fp, "# %*s:\t%s\n",        align * pad,  "MYOUT_FILEPATH",        args->myout_fileout );
   fprintf( fp, "# %*s:\t%s\n",        align * pad,  "MYDOMOUT_FILEPATH",     args->mydom_fileout );   
   fprintf( fp, "# ==============================================\n\n");
}

/*! FUNCTION:  ARGS_FindFiletype()
 *  SYNOPSIS:  Examines target and query, and finds the type of the files (by extension).
 */
FILETYPE 
ARGS_FindFiletype( STR    filename )
{
   for (int i = 0; i < NUM_FILETYPE_EXTS; i++) {
      char* ext_name = FILETYPE_EXTS[i].s;
      int   ext_type = FILETYPE_EXTS[i].i;
      if ( STRING_EndsWith( filename, ext_name, strlen(ext_name) ) == 0 ) {
         return ext_type;
      }
   }

   fprintf(stderr, "WARNING: '%s' filetype could not be found. Set to TYPE_NULL.\n", filename);
   
   return FILE_NULL;
}

/*! FUNCTION:  ARGS_Help_Info()
 *  SYNOPSIS:  Output help info.
 */
void 
ARGS_Help_Info()
{
   /* basic usage */
   // printf("Usage: ./fb-pruner <command> <target_hmm_file> <query_fasta_file>\n\n");

   /* command help */
   // printf("%-10s\t%-10s\t%s\n",
   //    "COMMANDS",
   //    "NUM_MAIN_ARGS",
   //    "DESC");
   // for (int i = 0; i < num_flag_cmds; i++) {
   //    printf("%-10s\t%-10d\t%-10s\t%s\n", 
   //       COMMAND_OPTS[i].long_flag,
   //       COMMAND_OPTS[i].num_args,
   //       DATATYPE_NAMES[COMMAND_OPTS[i].data_type],
   //       COMMAND_OPTS[i].desc );
   // }

   // /* option help */
   // printf("%-10s\t%-10s\t%-10s\t%s\n",
   //    "FLAG",
   //    "NUM_ARGS",
   //    "ARG_TYPE",
   //    "DESCRIPTION");
   // for (int i = 0; i < num_flag_cmds; i++) {
   //    printf("%-10s\t%-10d\t%-10s\t%s\n", 
   //       COMMAND_OPTS[i].long_flag,
   //       COMMAND_OPTS[i].num_args,
   //       DATATYPE_NAMES[COMMAND_OPTS[i].data_type],
   //       COMMAND_OPTS[i].desc );
   // }
   // printf("\n");
   // ERRORCHECK_exit(EXIT_SUCCESS);
}

/*! FUNCTION:  ARGS_Version_Info()
 *  SYNOPSIS:  Output version info.
 */
void 
ARGS_Version_Info()
{
   
}

/*! FUNCTION:  ARGS_OPTS_SetOptions()
 *  SYNOPSIS:  Output Default Options.
 */
STATUS_FLAG
ARGS_SetOptions(  ARGS*       args,
                  ARG_OPTS*   arg_opts )
{
   /* add all commandline options */
   PTR   arg_locs[]     = { &args->t_index_filein, &args->q_index_filein };
   INT   arg_dtypes[]   = { DATATYPE_STRING,    DATATYPE_STRING };

   {  ARG_OPTS_AddOption( arg_opts, "VERBOSE", 
         "Level of Output.",
         "Level of Output: [0] Minimal Output [1] Errors [2] Errors+Warnings [3] Maximal Output. ( Default: [1] )",
         "--verbose", "-v",
         1,
         (PTR[]){ &args->verbose_level },
         (INT[]){ DATATYPE_INT } ); 
   }
   {  ARG_OPTS_AddOption( arg_opts, "NUM_THREADS", 
         "Number of Threads.",
         "Level of Output: [0] Minimal Output [1] Errors [2] Errors+Warnings [3] Maximal Output. ( Default: [1] )",
         "--num-threads", NULL,
         1,
         (PTR[]){ &args->num_threads },
         (INT[]){ DATATYPE_INT } ); 
   }
   {  ARG_OPTS_AddOption( arg_opts, "ENFORCE WARNINGS", 
         "Should program terminate on warnings?",
         "Should program terminate on warnings? [0] Bypass warnings [1] Terminate on warnings.  ( Default: [1] )",
         "--enforce-warnings", "-w",
         1,
         (PTR[]){ &args->enforce_warnings },
         (INT[]){ DATATYPE_BOOL } ); }
   {  ARG_OPTS_AddOption( arg_opts, "INDEX FILE", 
         "Index file location.",
         "Index file location [query_index, target_index].  ( Default: NULL )",
         "--index", NULL,
         2,
         (PTR[]){ &args->t_index_filein, &args->q_index_filein },
         (INT[]){ DATATYPE_STRING, DATATYPE_STRING } ); }
   {  ARG_OPTS_AddOption( arg_opts, "USE LOCAL TOOLS", 
         "Should local or system tools be used?",
         "Should local or system tools be used? [0] Use system tools [1] Use local tools.  ( Default: [1] )",
         "--use-local-tools", NULL,
         1,
         (PTR[]){ &args->is_use_local_tools },
         (INT[]){ DATATYPE_BOOL } ); }
}

