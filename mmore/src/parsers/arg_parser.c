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

/* Parses Arguments from the command line */
void   ARGS_Parse( ARGS*   args,
                   int     argc, 
                   char*   argv[] )
{
   int   num_main_args     = 2; 
   char* flag              = NULL;
   /* required arguments */
   int   req_args    = 0;
   /* remaining arguments counter */
   int   args_rem    = argc - 1;
   /* current argument index */
   int   arg_cur     = 1;

   printf("NUM_ARGS: %d\n", argc);

   ARGS_SetDefaults(args);

   /* if no <command> argument given, run test case if in debug mode */
   if (argc <= 1) {
      printf("Usage: mmore <command> <target_hmm_file> <query_fasta_file>\n");
      printf("Hint: For more information, use '-h'");
      ERRORCHECK_exit(1);
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
         // printf("%s ?? %s\n", argv[1], PIPELINE_NAMES[i]);
         if ( STR_Compare( argv[1], PIPELINES[i].name ) == 0 ) {
            args->pipeline_mode = i;
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
   printf("NUM_MAIN_ARGS: %d\n", num_main_args);
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

   /* number of args based on pipeline */
   if ( num_main_args == 2 )
   {
      /* second arg is target */
      args->t_filepath = STR_Set( args->t_filepath, argv[2] );
      /* third arg is query */
      args->q_filepath = STR_Set( args->q_filepath, argv[3] );
   }
   if ( num_main_args == 5 )
   {
      /* second arg is target */
      args->t_filepath = STR_Set( args->t_filepath, argv[2] );
      /* third arg is query */
      args->q_filepath = STR_Set( args->q_filepath, argv[3] );
      /* fourth arg is mmseqs target (profile) */
      args->t_mmseqs_p_filepath = STR_Set( args->t_mmseqs_p_filepath, argv[4] );
      /* fourth arg is mmseqs target (sequence) */
      args->t_mmseqs_s_filepath = STR_Set( args->t_mmseqs_s_filepath, argv[5] );
      /* fourth arg is mmseqs query */
      args->q_mmseqs_filepath = STR_Set( args->q_mmseqs_filepath, argv[6] );
   }
   else if ( num_main_args == 1 )
   {
      /* second arg is target */
      args->t_filepath = STR_Set( args->t_filepath, argv[2] );
   }
   args_rem    -= num_main_args;
   arg_cur     += num_main_args;

   /* parse remaining flags and options */
   for (int i = arg_cur; i < argc; ++i)
   {
      /* if long flag */
      if ( STR_Compare_Prefix(argv[i], "--", 2) == 0 ) 
      {
         if ( STR_Compare( argv[i], (flag = "--help") ) == 0 ) {
            ARGS_Help_Info();
         }
         else if ( STR_Compare( argv[i], (flag = "--verbose") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->verbose_level = atoi(argv[i]);
               if ( args->verbose_level < 0 || args->verbose_level > NUM_VERBOSITY_MODES - 1 ) {
                  fprintf(stderr, "ERROR: Verbose level (%d) is outside acceptable range (%d,%d).\n", 
                     args->verbose_level, 0, NUM_VERBOSITY_MODES - 1 );
                  args->verbose_level = MAX( args->verbose_level, 0 );
                  args->verbose_level = MIN( args->verbose_level, NUM_VERBOSITY_MODES-1 );
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
         else if ( STR_Compare( argv[i], (flag = "--adjust-mmseqs-aln") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               if ( atoi(argv[i]) == 0 ) {
                  args->adjust_mmseqs_alns = false;
               } else if ( atoi(argv[i]) == 1 ) {
                  args->adjust_mmseqs_alns = true;
               }
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         /* === INPUT FILES === */
         else if ( STR_Compare( argv[i], (flag = "--index") ) == 0 ) {
            req_args = 2;
            if (i+req_args <= argc) {
               i++;
               args->t_indexpath = STR_Set( args->t_indexpath, argv[i] );
               i++;
               args->q_indexpath = STR_Set( args->q_indexpath, argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument=.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--local-tools") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               args->is_use_local_tools = atoi(argv[i]); 
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument=.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--guess-type") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               args->is_guess_filetype = atoi(argv[i]); 
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument=.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--mmore-ftype") ) == 0 ) {
            req_args = 5;
            if (i+req_args <= argc) {
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
               fprintf(stderr, "ERROR: %s flag requires (%d) argument=.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--mmore-main-ftype") ) == 0 ) {
            req_args = 2;
            if (i+req_args <= argc) {
               i++;
               args->t_filetype = atoi(argv[i]); 
               i++;
               args->q_filetype = atoi(argv[i]);
               /* since filetype supplied, turn off guesser */
               args->is_guess_filetype = false;
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument=.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--tmp") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               args->tmp_folderpath = STR_Set( args->tmp_folderpath, argv[i] );
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument=.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--mmseqs-m8") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               args->mmseqs_m8_filepath = STR_Set( args->mmseqs_m8_filepath, argv[i] );
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
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
         else if ( STR_Compare( argv[i], (flag = "--eval") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               args->mmseqs_evalue = atof(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--run-bias") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               args->is_compo_bias = atof(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--run-full") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               args->is_run_full = atof(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--run-domains") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               args->is_run_domains = atof(argv[i]);
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
               args->mmseqs_kmer = atof(argv[i]);
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
         /* === OUTPUT === */
         else if ( STR_Compare( argv[i], (flag = "--stderr") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               fclose(stdout);
               free(args->error_filepath);
               args->error_filepath = STR_Create(argv[i]);
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
               free(args->output_filepath);
               args->output_filepath = STR_Create(argv[i]);
               args->is_redirect_stdout = true;
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--domtblout") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               ERROR_free(args->tblout_filepath);
               args->tblout_filepath = STR_Create(argv[i]);
               args->is_tblout = true;
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--m8out") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               ERROR_free(args->m8out_filepath);
               args->m8out_filepath = STR_Create(argv[i]);
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
               ERROR_free(args->myout_filepath);
               args->myout_filepath = STR_Create(argv[i]);
               args->is_myout = true;
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--mydomout") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               ERROR_free(args->mydomout_filepath);
               args->mydomout_filepath = STR_Create(argv[i]);
               args->is_mydomout = true;
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               ERRORCHECK_exit(EXIT_FAILURE);
            }
         }
         else if ( STR_Compare( argv[i], (flag = "--mytimeout") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               ERROR_free(args->mytimeout_filepath);
               args->mytimeout_filepath = STR_Create(argv[i]);
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
               ERROR_free(args->mythreshout_filepath);
               args->mythreshout_filepath = STR_Create(argv[i]);
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
               args->output_filepath = STR_Create(argv[i]);
               args->is_redirect_stdout = true;
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
         else if ( STR_Compare( argv[i], (flag = "--customout") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               ERROR_free(args->tblout_filepath);
               args->tblout_filepath = STR_Create(argv[i]);
               args->is_tblout = true;
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

   // if ( args->is_guess_filetype == true ) {
   //    args->t_filetype           = ARGS_Find_FileType( args->t_filepath );
   //    args->q_filetype           = ARGS_Find_FileType( args->q_filepath );
   //    args->t_mmseqs_p_filetype    = ARGS_Find_FileType( args->t_mmseqs_p_filepath );
   //    args->q_mmseqs_filetype    = ARGS_Find_FileType( args->q_mmseqs_filepath );
   // }
}

/* SET DEFAULT ARGUMENTS (generic) */
void  ARGS_SetDefaults( ARGS* args )
{
   /* --- PIPELINE OPTIONS --- */
   args->pipeline_mode           = PIPELINE_TEST;
   args->verbose_level           = VERBOSE_LOW;
   args->search_mode             = MODE_UNILOCAL;
   args->qt_search_space         = SELECT_ALL_V_ALL;
   args->tmp_folderpath          = NULL;
   args->tmp_remove              = false;

   /* --- TASK OPTIONS --- */
   args->is_compo_bias           = true;
   args->is_run_pruned           = true;
   args->is_run_full             = false;
   args->is_run_domains          = true;

   /* --- DEBUG OPTIONS --- */
   args->is_use_local_tools      = false;
   args->is_recycle_mx           = false;
   args->is_debug                = true;
   args->dbg_folderpath          = STR_Create("test_output/");

   /* --- INPUT --- */
   /* filepath */
   args->t_filepath              = NULL;
   args->q_filepath              = NULL;
   args->t_mmseqs_p_filepath     = NULL;
   args->t_mmseqs_p_filepath     = NULL;
   args->q_mmseqs_filepath       = NULL;
   /* filetype */
   args->is_guess_filetype       = true;
   args->t_filetype              = FILE_HMM;
   args->q_filetype              = FILE_FASTA;
   args->t_mmseqs_p_filetype     = FILE_HHM;
   args->t_mmseqs_p_filetype     = FILE_HHM;
   args->t_mmseqs_p_filetype     = FILE_HHM;
   /* indexes */
   args->t_indexpath             = NULL;
   args->q_indexpath             = NULL;

   /* --- INTERRIM OUTPUT --- */
   args->mmseqs_m8_filepath     = STR_Create("mmseqs.results.m8out");

   /* --- OUTPUT --- */
   args->is_redirect_stdout      = false;
   args->output_filepath         = STR_Create("results.stdout");
   args->is_tblout               = false;
   args->tblout_filepath         = STR_Create("results.tblout");
   args->is_m8out                = true;
   args->m8out_filepath          = STR_Create("results.m8");
   args->is_myout                = false;
   args->myout_filepath          = STR_Create("results.myout");
   args->is_mydomout             = false;
   args->mydomout_filepath       = STR_Create("results.mydomout");
   args->is_mytimeout            = false;
   args->mytimeout_filepath      = STR_Create("results.mytimeout");
   args->is_mythreshout          = false;
   args->mythreshout_filepath    = STR_Create("results.mythreshout");

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
   args->mmseqs_kmer             = 7;
   args->mmseqs_prefilter        = 80;
   args->mmseqs_ungapped_vit     = 15;
   args->mmseqs_evalue           = 1000.0;
   args->mmseqs_pvalue           = 1e-3f;

   /* --- VITERBI / FWDBACK THRESHOLDS --- */
   args->is_run_filter           = false;
   args->threshold_vit           = 1e-3f;
   args->threshold_cloud         = 1e-5f;
   args->threshold_bound_fwd     = 1e-5f;
   args->threshold_fwd           = 1e-5f;
   args->threshold_mmore         = 1e-5f;
}


/* sends ARGS data to FILE POINTER */
void ARGS_Dump( ARGS*    args,
                FILE*    fp )
{
   int      pad               = 30;
   bool     align             = 1;     /* -1 for right alignment, 1 for left alignment */

   fprintf( fp, "=== ARGS =====================\n");
   /* --- PIPELINE --- */
   fprintf( fp, "%*s:\t%s == %d\n",    align * pad,  "PIPELINE",        PIPELINES[args->pipeline_mode].name,   args->pipeline_mode );
   fprintf( fp, "%*s:\t%s == %d\n",    align * pad,  "VERBOSITY_MODE",  VERBOSITY_NAMES[args->verbose_level],  args->verbose_level );
   fprintf( fp, "%*s:\t%s\n",          align * pad,  "SEARCH_MODE",     MODE_NAMES[args->search_mode] );
   fprintf( fp, "%*s:\t%s\n",          align * pad,  "COMPO_BIAS",      ( args->is_compo_bias ? "True" : "False" ) );
   fprintf( fp, "\n" );
   /* --- INPUT --- */
   fprintf( fp, "%*s:\t%s\n",          align * pad,  "TARGET_FILEPATH", args->t_filepath );
   fprintf( fp, "%*s:\t%s\n",          align * pad,  "TARGET_FILETYPE", FILE_TYPE_NAMES[args->t_filetype] );
   fprintf( fp, "%*s:\t%s\n",          align * pad,  "QUERY_FILEPATH",  args->q_filepath );
   fprintf( fp, "%*s:\t%s\n",          align * pad,  "QUERY_FILETYPE",  FILE_TYPE_NAMES[args->q_filetype] );
   fprintf( fp, "%*s:\t%s\n",          align * pad,  "TARGET_MMSEQS_FILEPATH",  args->t_mmseqs_p_filepath );
   fprintf( fp, "%*s:\t%s\n",          align * pad,  "TARGET_MMSEQS_FILETYPE",  FILE_TYPE_NAMES[args->t_mmseqs_p_filetype] );
   fprintf( fp, "%*s:\t%s\n",          align * pad,  "T_INDEX_PATH",    args->t_indexpath );
   fprintf( fp, "%*s:\t%s\n",          align * pad,  "Q_INDEX_PATH",    args->q_indexpath );
   fprintf( fp, "%*s:\t%s\n",          align * pad,  "MMSEQS_RESULTS",  args->mmseqs_m8_filepath );
   fprintf( fp, "%*s:\t%s\n",          align * pad,  "TMP_FOLDERPATH",  args->tmp_folderpath ); 
   fprintf( fp, "\n" );
   /* --- MMSEQS --- */
   fprintf( fp, "%*s:\t%d\n",          align * pad,  "MMSEQS_KMER",     args->mmseqs_kmer );
   fprintf( fp, "%*s:\t%d\n",          align * pad,  "MMSEQS_KSCORE",   args->mmseqs_prefilter );
   fprintf( fp, "%*s:\t%d\n",          align * pad,  "MMSEQS_UNGAPPED", args->mmseqs_ungapped_vit );
   fprintf( fp, "%*s:\t%.3g\n",        align * pad,  "MMSEQS_EVALUE",   args->mmseqs_evalue );
   fprintf( fp, "%*s:\t%.3g\n",        align * pad,  "MMSEQS_PVALUE",   args->mmseqs_pvalue ); 
   /* --- MMORE / FB-PRUNER --- */ 
   fprintf( fp, "%*s:\t%.3g\n",        align * pad,  "ALPHA",           args->alpha );
   fprintf( fp, "%*s:\t%.3g\n",        align * pad,  "BETA",            args->beta );
   fprintf( fp, "%*s:\t%d\n",          align * pad,  "GAMMA",           args->gamma );
   fprintf( fp, "%*s:\t%.3g\n",        align * pad,  "E_VALUE",         args->mmore_evalue );
   fprintf( fp, "%*s:\t%.3g\n",        align * pad,  "P_VALUE",         args->mmore_pvalue );
   fprintf( fp, "\n" );
   /* --- OUTPUT --- */
                           fprintf( fp, "%*s:\t%s\n",          align * pad,  "OUTPUT_FILEPATH", args->output_filepath );
   if (args->is_tblout)    fprintf( fp, "%*s:\t%s\n",          align * pad,  "TBLOUT_FILEPATH", args->tblout_filepath );
   if (args->is_m8out)     fprintf( fp, "%*s:\t%s\n",          align * pad,  "M8OUT_FILEPATH", args->m8out_filepath );
   if (args->is_myout)     fprintf( fp, "%*s:\t%s\n",          align * pad,  "MYOUT_FILEPATH", args->myout_filepath );
   if (args->is_mydomout)  fprintf( fp, "%*s:\t%s\n",          align * pad,  "MYDOMOUT_FILEPATH", args->mydomout_filepath );
   if (args->is_customout) fprintf( fp, "%*s:\t%s\n",          align * pad,  "CUSTOMOUT_FILEPATH", args->customout_filepath );
   
   fprintf( fp, "=============================\n\n");
}

/* examines target and query, and finds the type of the files */
FILE_TYPE 
ARGS_Find_FileType( char* filename )
{
   for (int i = 0; i < NUM_FILE_EXTS; i++) {
      char* ext = FILE_TYPE_EXTS[i];
      if ( STRING_EndsWith( filename, ext, strlen(ext) ) == 0 ) {
         return FILE_TYPE_MAP[i];
      }
   }

   fprintf(stderr, "WARNING: '%s' filetype could not be found. Set to TYPE_NULL.\n", filename);
   
   return FILE_NULL;
}

/* output help info */
void 
ARGS_Help_Info()
{
   /* basic usage */
   printf("Usage: ./fb-pruner <command> <target_hmm_file> <query_fasta_file>\n\n");

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

   /* option help */
   printf("%-10s\t%-10s\t%-10s\t%s\n",
      "FLAG",
      "NUM_ARGS",
      "ARG_TYPE",
      "DESCRIPTION");
   for (int i = 0; i < num_flag_cmds; i++) {
      printf("%-10s\t%-10d\t%-10s\t%s\n", 
         COMMAND_OPTS[i].long_flag,
         COMMAND_OPTS[i].num_args,
         DATATYPE_NAMES[COMMAND_OPTS[i].data_type],
         COMMAND_OPTS[i].desc );
   }
   printf("\n");
   ERRORCHECK_exit(EXIT_SUCCESS);
}

/* output version info */
void ARGS_Version_Info()
{
   
}

/* initialize default arg options */
STATUS_FLAG
ARG_OPTS_Default_Opts( ARGS* args )
{
   /* command line flags and options */
   int   num_flag_cmds = 11;
   // ARG_OPT COMMAND_OPTS[] = {
   //    /* name | num_args | data_type | arg_loc | arg_bool | long_flag | short_flag | desc */
   //    /* output formats */
   //    {  "OUTFILE",           1,    DATATYPE_INT,        NULL,    "--output",          "-o",    "Result output file destination [test_output/results.tsv]."  },
   //    {  "INFILES",           2,    DATATYPE_STRING,     NULL,    "--input",           "-i",    "Input files: {target,query} [test cases]."  },
   //    /* general options */
   //    {  "INDEX",             2,    DATATYPE_STRING,     NULL,    "--index",           "-x",    "Index files: {target,query} [builds on fly]."  },
   //    /* mmseqs options */
   //    {  "MMSEQS_TMP",        1,    DATATYPE_STRING,     NULL,    "--mmseqs-tmp",      NULL,    "MMseqs temp folder [null]."  },
   //    {  "MMSEQS_INPUT",      1,    DATATYPE_STRING,     NULL,    "--mmseqs-input",    NULL,    "MMseqs results file input [null]."  },
   //    {  "MMSEQS_LOOKUP",     2,    DATATYPE_STRING,     NULL,    "--mmseqs-lookup",   NULL,    "MMseqs lookup files: {target,query} [null]."  },
   //    /* mmore main options */
   //    {  "ALPHA",             1,    DATATYPE_FLOAT,      NULL,    "--alpha",           "-a",    "MMORE X-drop per antidiagonal pruning ratio [20.0]." },
   //    {  "BETA",              1,    DATATYPE_INT,        NULL,    "--beta",            "-b",    "MMore X-drop global " },
   //    {  "BETA",              1,    DATATYPE_INT,        NULL,    "--beta",            "-b",    "Number of passes of cloud search before pruning [5]." },
      
   //    {  "WINDOW",            4,    DATATYPE_INT,        NULL,    "--window",          "-w",    "Examine substring of query and target."  },
   //    {  "Q_RANGE",           2,    DATATYPE_INT,        NULL,    "--qrange",          NULL,    "Give range of ids in query file index to search [-1,-1]."  },
   //    {  "T_RANGE",           2,    DATATYPE_INT,        NULL,    "--trange",          NULL,    "Give range of ids in target file index to search [-1,-1]."  },
   // };
}