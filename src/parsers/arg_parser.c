/*******************************************************************************
 *     - FILE:   arg_parser.h
 *  - DESC:    Parses command line arguments.
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

void ARGS_MainArg_Parser(ARGS* args,
                         COMMANDLINE* cmd);

void ARGS_Parse(ARGS* args,
                int argc,
                char* argv[],
                COMMANDLINE* cmd,
                ARG_OPTS* arg_opts) {
  int num_main_args = 2;
  char* flag = NULL;
  /* required arguments */
  int req_args = 0;
  /* remaining arguments counter */
  int args_rem = argc - 1;
  /* current argument index */
  int arg_cur = 1;

  DBG_PRINTF(stdout, "# === NUM_ARGS: %d\n", argc);
  ARGS_SetDefaults(args);
  // ARGS_SetOptions( args, arg_opts );

  /* if no <command> argument given, run test case if in debug mode */
  if (argc <= 1) {
    printf("Usage: mmoreseqs <command> <main_args...> <opts...>\n");
    printf("Hint: For more information, use '-h'");
    ERRORCHECK_exit(EXIT_FAILURE);
  }

  /* check for help flag */
  if (STR_Equals(argv[1], "-h") || STR_Equals(argv[1], "--help")) {
    ARGS_Help_Info();
  }

  /* check for version flag */
  if (STR_Equals(argv[1], "--version")) {
    ARGS_Version_Info();
  }

  /* first argument is command pipeline (id_number or name) */
  bool found_pipeline = false;
  for (int i = 0; i < NUM_PIPELINE_MODES; i++) {
    PIPELINE* pipeline = &PIPELINES[i];
    if (STR_Compare(argv[1], pipeline->name) == 0) {
      args->pipeline_mode = i;
      args->pipeline_name = STR_Create(pipeline->name);
      found_pipeline = true;
      break;
    }
  }

  /* check that valid pipeline mode was entered */
  if (found_pipeline == false || (args->pipeline_mode < 0) || (args->pipeline_mode >= NUM_PIPELINE_MODES)) {
    fprintf(stderr, "ERROR: Invalid pipeline/command was given: (%s, %d).\n", argv[1], args->pipeline_mode);
    fprintf(stderr, "VALID PIPELINE/COMMANDS OPTS: [ ");
    for (int i = 0; i < NUM_PIPELINE_MODES; i++)
      fprintf(stderr, "%s, ", PIPELINES[i].name);
    fprintf(stderr, "]\n");
    ERRORCHECK_exit(EXIT_FAILURE);
  }
  args_rem -= 1;
  arg_cur += 1;

  /* set number of main arguments based on given pipeline */
  num_main_args = PIPELINES[args->pipeline_mode].num_main_args;
  /* check proper number of main args remain */
  if (args_rem < num_main_args) {
    fprintf(stderr, "ERROR: Improper number of main arguments. [required: %d/%d]\n", args_rem, num_main_args);
    ARGS_Command_Help_Info(args);
#if DEBUG
    printf("Using DEFAULT arguments...\n\n");
    return;
#else
    ERRORCHECK_exit(EXIT_FAILURE);
#endif
  }

  /* check if proper number of args */
  for (int i = 0; i < num_main_args; i++) {
    if (STR_StartsWith(argv[i + 1], "--") == 0) {
      fprintf(stderr, "ERROR: Improper number of main arguments. [required: %d/%d]\n", i, num_main_args);
      ARGS_Command_Help_Info(args);
      ERRORCHECK_exit(EXIT_FAILURE);
    }
  }

  /* parse main commands */
  if (STR_Equals(args->pipeline_name, "search")) {
    args->t_mmore_filein = STR_Set(args->t_mmore_filein, argv[2]);
    args->q_mmore_filein = STR_Set(args->q_mmore_filein, argv[3]);
    args->t_mmseqs_p_filein = STR_Set(args->t_mmseqs_p_filein, argv[4]);
    args->t_mmseqs_s_filein = STR_Set(args->t_mmseqs_s_filein, argv[5]);
    args->q_mmseqs_filein = STR_Set(args->q_mmseqs_filein, argv[6]);

    args->t_mmore_p_filetype = FILE_HMM;
    args->q_filetype = FILE_FASTA;
    args->t_mmseqs_p_filetype = FILE_MMDB_P;
    args->t_mmseqs_s_filetype = FILE_MMDB_S;
    args->q_mmseqs_filetype = FILE_MMDB_S;
  }
  elif (STR_Equals(args->pipeline_name, "mmore-search")) {
    args->t_filein = STR_Set(args->t_filein, argv[2]);
    args->q_filein = STR_Set(args->q_filein, argv[3]);
    args->mmseqs_m8_filein = STR_Set(args->mmseqs_m8_filein, argv[4]);

    args->t_filetype = FILE_HMM;
    args->q_filetype = FILE_FASTA;
  }
  elif (STR_Equals(args->pipeline_name, "mmseqs-search")) {
    args->t_mmseqs_p_filein = STR_Set(args->t_mmseqs_p_filein, argv[2]);
    args->t_mmseqs_s_filein = STR_Set(args->t_mmseqs_s_filein, argv[3]);
    args->q_mmseqs_filein = STR_Set(args->q_mmseqs_filein, argv[4]);

    args->t_mmseqs_p_filetype = FILE_MMDB_P;
    args->t_mmseqs_s_filetype = FILE_MMDB_S;
    args->q_mmseqs_filetype = FILE_MMDB_S;
  }
  elif (STR_Equals(args->pipeline_name, "prep")) {
    args->target_prep = STR_Set(args->target_prep, argv[2]);
    args->query_prep = STR_Set(args->query_prep, argv[3]);
    args->tmp_folderpath = STR_Set(args->tmp_folderpath, argv[4]);
    args->prep_folderpath = STR_Set(args->prep_folderpath, argv[4]);

    args->target_prep_type = FILE_MSA;
    args->query_prep_type = FILE_FASTA;
  }
  elif (STR_Equals(args->pipeline_name, "prep-search")) {
    args->prep_folderpath = STR_Set(args->prep_folderpath, argv[2]);
    args->tmp_folderpath = STR_Set(args->tmp_folderpath, argv[2]);
  }
  elif (STR_Equals(args->pipeline_name, "easy-search")) {
    args->target_prep = STR_Set(args->target_prep, argv[2]);
    args->query_prep = STR_Set(args->query_prep, argv[3]);
    args->tmp_folderpath = STR_Set(args->tmp_folderpath, argv[4]);
    args->prep_folderpath = STR_Set(args->prep_folderpath, argv[4]);

    args->target_prep_type = FILE_MSA;
    args->query_prep_type = FILE_FASTA;
  }
  elif (STR_Equals(args->pipeline_name, "index")) {
    args->t_index_filein = STR_Set(args->t_index_filein, argv[2]);
    // args->t_indexout     = STR_Set( args->t_indexout,     argv[3] );
  }
  else {
    fprintf(stderr, "ERROR: Pipeline option '%s' is currently not supported.\n", args->pipeline_name);
    ERRORCHECK_exit(EXIT_FAILURE);
  }
  args_rem -= num_main_args;
  arg_cur += num_main_args;

  /* parse remaining flags and options */
  for (int i = arg_cur; i < argc; ++i) {
    /* if long flag */
    if (STR_ComparePrefix(argv[i], "--", 2) == 0) {
      /* === HELP OPTIONS === */
      if (STR_Equals(argv[i], (flag = "--help"))) {
        ARGS_Command_Help_Info(args);
        ERRORCHECK_exit(EXIT_SUCCESS);
      }
      elif (STR_Equals(argv[i], (flag = "--info"))) {
        ARGS_Version_Info();
      }
      elif (STR_Equals(argv[i], (flag = "--version"))) {
        ARGS_Version_Info();
      }
      /* === DEBUG OPTIONS (should only affect debug builds) === */
      elif (STR_Equals(argv[i], (flag = "--dbg"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          debugger->is_debugging = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--dbg-viz"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          debugger->is_viz = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === PIPELINE OPTIONS === */
      elif (STR_Equals(argv[i], (flag = "--verbose"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->verbose_level = atoi(argv[i]);
          if (args->verbose_level < 0 || args->verbose_level > NUM_VERBOSITY_MODES - 1) {
            fprintf(stderr, "ERROR: Verbose Level (%d) is outside acceptable range (%d,%d).\n",
                    args->verbose_level, 0, NUM_VERBOSITY_MODES - 1);
            args->verbose_level = MAX(args->verbose_level, 0);
            args->verbose_level = MIN(args->verbose_level, NUM_VERBOSITY_MODES - 1);
            fprintf(stderr, "WARNING: Verbose Level set to: %d.\n", args->verbose_level);
          }
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--num-threads"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->num_threads = atoi(argv[i]);
          args->num_threads = atoi(argv[i]);
          if (args->num_threads < 1) {
            fprintf(stderr, "ERROR: Number of Threads (%d) is outside acceptable range (%d,%d).\n",
                    args->num_threads, 1, 16);
            args->num_threads = MAX(args->num_threads, 1);
            args->num_threads = MIN(args->num_threads, 16);
            fprintf(stderr, "WARNING: Number of Threads set to: %d.\n", args->num_threads);
          }
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--enforce-errors"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          if (atoi(argv[i])) {
            args->enforce_warnings = false;
          }
          elif (atoi(argv[i]) == 1) {
            args->enforce_warnings = true;
          }
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--enforce-warnings"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          if (atoi(argv[i])) {
            args->enforce_warnings = false;
          }
          elif (atoi(argv[i]) == 1) {
            args->enforce_warnings = true;
          }
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--search-type"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          if (STR_Equals(argv[i], "P2S") == true) {
            args->search_name = STR_Set(args->search_name, argv[i]);
            args->target_prep_type = FILE_MSA;
            args->query_prep_type = FILE_FASTA;
          }
          elif (STR_Equals(argv[i], "S2S") == true) {
            args->search_name = STR_Set(args->search_name, argv[i]);
            args->target_prep_type = FILE_FASTA;
            args->query_prep_type = FILE_FASTA;
          }
          else {
            fprintf(stderr, "ERROR: '%s' is not a valid argument for --search-type.\n", argv[i]);
            ERRORCHECK_exit(EXIT_FAILURE);
          }
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === TASKS === */
      elif (STR_Equals(argv[i], (flag = "--run-mmseqs"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_mmseqs = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-mmseqs-pref"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_mmseqs_pref = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-mmseqs-align"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_mmseqs_align = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-convert"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_convert = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-mmore"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_mmore = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--use-pvals"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->use_pvals = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === INPUT PROGRAMS === */
      elif (STR_Equals(argv[i], (flag = "--program-mmseqs"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmseqs_program = STR_Set(args->mmseqs_program, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--program-hmmer"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->hmmer_program = STR_Set(args->hmmer_program, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--program-mmoreseqs"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmoreseqs_program = STR_Set(args->mmoreseqs_program, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--script-dir"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmoreseqs_scripts = STR_Set(args->mmoreseqs_scripts, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === INPUT FILES === */
      elif (STR_Equals(argv[i], (flag = "--index"))) {
        req_args = 2;
        if (i + req_args < argc) {
          i++;
          args->t_index_filein = STR_Set(args->t_index_filein, argv[i]);
          i++;
          args->q_index_filein = STR_Set(args->q_index_filein, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument=.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--local-tools"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_use_local_tools = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument=.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--guess-ftype"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_guess_filetype = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument=.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmoreseqs-ftype"))) {
        req_args = 5;
        if (i + req_args < argc) {
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
      elif (STR_Equals(argv[i], (flag = "--mmoreseqs-main-ftype"))) {
        req_args = 2;
        if (i + req_args < argc) {
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
      elif (STR_Equals(argv[i], (flag = "--tmp"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->tmp_folderpath = STR_Set(args->tmp_folderpath, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--prep"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->tmp_folderpath = STR_Set(args->prep_folderpath, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmseqs-m8"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmseqs_m8_filein = STR_Set(args->mmseqs_m8_filein, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === PREP FILES === */
      elif (STR_Equals(argv[i], (flag = "--prep-link-target-mmore"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->link_target_mmore_prep = STR_Set(args->link_target_mmore_prep, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--prep-link-query-mmore"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->link_query_mmore_prep = STR_Set(args->link_query_mmore_prep, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--prep-link-target-mmseqs"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->link_target_mmseqs_prep = STR_Set(args->link_target_mmseqs_prep, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--prep-link-query-mmseqs"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->link_query_mmseqs_prep = STR_Set(args->link_query_mmseqs_prep, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === MMORE METADATA === */
      elif (STR_Equals(argv[i], (flag = "--dbsizes"))) {
        req_args = 2;
        if (i + req_args < argc) {
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
      elif (STR_Equals(argv[i], (flag = "--alpha"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->alpha = atof(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--beta"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->beta = atof(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--gamma"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->gamma = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--hard-limit"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->hard_limit = atof(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* ==== MMORE OPTIONS === */
      elif (STR_Equals(argv[i], (flag = "--run-prep"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->tmp_folderpath = STR_Set(args->tmp_folderpath, argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-bias"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_bias = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-full"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_full = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-domains"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_domains = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-vit-mmore"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_vit_mmore = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-mmseqsaln"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_mmseqsaln = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-vitaln"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_vit = atoi(argv[i]);
          args->is_run_vitaln = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-vit"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_vit = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-postaln"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_optacc = atoi(argv[i]);
          args->is_run_postaln = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-post"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_optacc = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === MMORE FILTERS === */
      elif (STR_Equals(argv[i], (flag = "--run-filter"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_filter = atoi(argv[i]);
          args->is_run_viterbi_filter = atoi(argv[i]);
          args->is_run_cloud_filter = atoi(argv[i]);
          args->is_run_boundfwd_filter = atoi(argv[i]);
          args->is_run_fwdback_filter = atoi(argv[i]);
          args->is_run_report_filter = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-vit-filter"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_viterbi_filter = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-cld-filter"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_cloud_filter = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--run-fwd-filter"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->is_run_boundfwd_filter = atoi(argv[i]);
          args->is_run_fwdback_filter = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--vit-filter"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->threshold_vit = atof(argv[i]);
          // args->is_run_viterbi_filter   = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--cld-filter"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->threshold_cloud = atof(argv[i]);
          // args->is_run_cloud_filter  = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--fwd-filter"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->threshold_boundfwd = atof(argv[i]);
          args->threshold_fwd = atof(argv[i]);
          // args->is_run_boundfwd_filter  = true;
          // args->is_run_fwdback_filter   = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--eval"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->threshold_report_eval = atof(argv[i]);
          args->is_run_report_filter = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === MMSEQS PARAMETERS === */
      elif (STR_Equals(argv[i], (flag = "--mmseqs-split"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmseqs_split = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmseqs-kmer"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmseqs_kmer = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmseqs-kscore"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmseqs_kscore = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmseqs-sens"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmseqs_sensitivity = atof(argv[i]);
          args->is_run_mmseqs_sens = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmseqs-ungapped-vit"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmseqs_ungapped_vit = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmseqs-eval"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmseqs_evalue = atof(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmseqs-pval"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmseqs_pvalue = atof(argv[i]);
          args->mmseqs_p2s_pvalue = atof(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmseqs-hits-per-search"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmseqs_maxhits = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmseqs-altalis"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->mmseqs_altalis = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === MMSEQS METADATA === */
      elif (STR_Equals(argv[i], (flag = "--mmseqs-times"))) {
        req_args = 5;
        if (i + req_args < argc) {
          i++;
          args->prep_time = atof(argv[i]);
          i++;
          args->mmseqs_time = atof(argv[i]);
          i++;
          args->mmseqs_prefilter_time = atof(argv[i]);
          i++;
          args->mmseqs_p2s_time = atof(argv[i]);
          i++;
          args->mmseqs_s2s_time = atof(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mmseqs-dbsizes"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->prefilter_dbsize = atoi(argv[i]);
          i++;
          args->mmseqs_dbsize = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === SEARCH/RANGE OPTIONS === */
      elif (STR_Equals(argv[i], (flag = "--range"))) {
        req_args = 2;
        if (i + req_args < argc) {
          i++;
          args->list_range.beg = atoi(argv[i]);
          i++;
          args->list_range.end = atoi(argv[i]);
        } else {
          fprintf(stderr, "ERROR: '%s' flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--search-mode"))) {
        req_args = 1;
        if (i + req_args < argc) {
          i++;
          args->search_mode = atoi(argv[i]);
          if (args->search_mode < 0 || args->search_mode > NUM_SELECT_SEARCHES) {
            fprintf(stderr, "ERROR: Verbose level (%d) is outside acceptable range (%d,%d).\n",
                    args->search_mode, 0, NUM_SELECT_SEARCHES);
            args->search_mode = MAX(args->search_mode, 0);
            args->search_mode = MIN(args->search_mode, NUM_SELECT_SEARCHES - 1);
          }
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === INTERRIM OUTPUT === */
      elif (STR_Equals(argv[i], (flag = "--mmseqs-m8out"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          free(args->mmseqs_m8out);
          args->mmseqs_m8out = STR_Create(argv[i]);
          args->is_mmseqs_m8out = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      /* === OUTPUT === */
      elif (STR_Equals(argv[i], (flag = "--stderr"))) {
        req_args = 1;
        if (i + req_args <= argc) {
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
      elif (STR_Equals(argv[i], (flag = "--stdout"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          free(args->stdout_fileout);
          args->stdout_fileout = STR_Create(argv[i]);
          args->is_redirect_stdout = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--allout"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          STR_Destroy(args->hmmerout_fileout);
          args->hmmerout_fileout = STR_Concat(argv[i], ".hmmerout");
          STR_Destroy(args->m8out_fileout);
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
      elif (STR_Equals(argv[i], (flag = "--domtblout"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          ERROR_free(args->mydom_fileout);
          args->mydom_fileout = STR_Create(argv[i]);
          args->is_mydom = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--m8out"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          ERROR_free(args->m8out_fileout);
          args->m8out_fileout = STR_Create(argv[i]);
          args->is_m8out = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--myout"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          ERROR_free(args->myout_fileout);
          args->myout_fileout = STR_Create(argv[i]);
          args->is_myout = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mydomtblout"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          ERROR_free(args->mydom_fileout);
          args->mydom_fileout = STR_Create(argv[i]);
          args->is_mydom = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mytimeout"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          ERROR_free(args->mytime_fileout);
          args->mytime_fileout = STR_Create(argv[i]);
          args->is_mytimeout = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--mythreshout"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          ERROR_free(args->mythresh_fileout);
          args->mythresh_fileout = STR_Create(argv[i]);
          args->is_mythreshout = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--customout"))) {
        req_args = 1;
        if (i + req_args <= argc) {
          i++;
          ERROR_free(args->hmmerout_fileout);
          args->hmmerout_fileout = STR_Create(argv[i]);
          args->is_hmmerout = true;
        } else {
          fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
          ERRORCHECK_exit(EXIT_FAILURE);
        }
      }
      elif (STR_Equals(argv[i], (flag = "--debugout"))) {
        req_args = 1;
        if (i + req_args <= argc) {
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
    } else {
      fprintf(stderr, "ERROR: '%s' is not associated with an argument or flag\n", argv[i]);
      ERRORCHECK_exit(EXIT_FAILURE);
    }
  }
}

void ARGS_MainArg_Parser(ARGS* args,
                         COMMANDLINE* cmd) {
}

bool ARGS_OptArg_Parser(ARGS* args,
                        char* argv[],
                        int argc,
                        void* arg_opts[],
                        char* in_flag,
                        char* flag,
                        int arg_req,
                        int arg_cur) {
  if (STR_Equals(in_flag, flag)) {
    if (arg_cur + arg_req < argc) {
      for (int i = 0; i < arg_req; i++) {
        arg_cur++;
        arg_opts[i] = STR_Set(arg_opts[i], argv[arg_cur]);
      }
    } else {
      fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, arg_req);
      ERRORCHECK_exit(EXIT_FAILURE);
    }
    return true;
  }
  return false;
}

void ARGS_SetDefaults(ARGS* args) {
  /* --- COMMANDLINE --- */

  /* --- PROGRAM --- */
  args->mmoreseqs_program = STR_Create("mmoreseqs");
  args->mmseqs_program = STR_Create("mmseqs");
  args->hmmer_program = STR_Create("hmmbuild");
  args->mmoreseqs_scripts = STR_Create(SCRIPT_DIR);
  /* --- PIPELINE MODE --- */
  args->pipeline_mode = PIPELINE_TEST;
  args->pipeline_name = NULL;
  args->verbose_level = VERBOSE_LOW;
  args->num_threads = 1;
  args->search_mode = MODE_UNILOCAL;
  args->qt_search_space = SELECT_ALL_V_ALL;
  args->tmp_folderpath = NULL;
  args->tmp_remove = false;

  /* --- PREPARATION OPTIONS --- */
  /* files/folders */
  args->is_run_prep = true;
  args->is_prep_copy = true;
  args->prep_folderpath = NULL;
  args->target_prep = NULL;
  args->query_prep = NULL;
  /* types */
  args->target_prep_type = FILE_MSA;
  args->query_prep_type = FILE_FASTA;

  /* --- PIPELINE OPTIONS --- */
  args->search_name = STR_Create("P2S");
  args->is_run_bias = true;
  args->is_run_pruned = true;
  args->is_run_full = false; /* DEBUG */
  args->is_run_domains = true;
  args->is_run_mmseqsaln = false;
  args->is_run_vit_mmore = false; /* DEBUG */
  args->is_run_vit = false;       /* DEBUG */
  args->is_run_vitaln = false;    /* DEBUG */
  args->is_run_optacc = false;    /* DEBUG */
  args->is_run_post = false;
  args->is_run_postaln = false;
  args->is_run_mmseqs = true;
  args->is_run_mmseqs_pref = true;
  args->is_run_mmseqs_align = true;
  args->is_run_mmseqs_ungapped = true;
  args->is_run_mmseqs_sens = false;
  args->is_run_mmore = true;
  args->is_run_convert = true;
  args->use_pvals = true;

  /* --- DEBUG OPTIONS --- */
  args->is_use_local_tools = false;
  args->is_recycle_mx = false;
  args->is_debug = true;
  args->dbg_folderpath = STR_Create(DEBUG_FOLDER "/");

  /* --- INPUT --- */
  /* filepath */
  args->t_filein = NULL;
  args->q_filein = NULL;
  args->t_mmore_filein = NULL;
  args->q_mmore_filein = NULL;
  args->t_mmseqs_p_filein = NULL;
  args->t_mmseqs_s_filein = NULL;
  args->q_mmseqs_filein = NULL;
  /* filetype */
  args->is_guess_filetype = true;
  args->t_filetype = FILE_HMM;
  args->q_filetype = FILE_FASTA;
  args->t_mmore_filetype = FILE_HMM;
  args->q_mmore_filetype = FILE_FASTA;
  args->t_mmseqs_p_filetype = FILE_MMDB_P;
  args->t_mmseqs_s_filetype = FILE_MMDB_S;
  args->q_mmseqs_filetype = FILE_MMDB_S;
  /* indexes */
  args->t_index_filein = NULL;
  args->q_index_filein = NULL;
  /* db sizes */
  args->q_dbsize = -1;
  args->t_dbsize = -1;

  /* --- INTERRIM OUTPUT --- */
  args->mmseqs_m8_filein = STR_Create("mmseqs.results.m8out");
  args->is_mmseqs_m8out = false;
  args->mmseqs_m8out = STR_Create("mmseqs.results.m8out");

  /* --- OUTPUT --- */
  /* default outputs */
  args->is_redirect_stdout = false;
  args->stdout_fileout = STR_Create("mmore.results.stdout");
  args->is_redirect_stderr = false;
  args->stderr_fileout = STR_Create("mmore.results.stderr");
  /* special outputs */
  args->is_allout = false;
  args->allout_fileout = STR_Create("mmore.results");
  args->is_hmmerout = false;
  args->hmmerout_fileout = STR_Create("mmore.results.tblout");
  args->is_m8out = false;
  args->m8out_fileout = STR_Create("mmore.results.m8");
  args->is_myout = false;
  args->myout_fileout = STR_Create("mmore.results.myout");
  args->is_mydom = false;
  args->mydom_fileout = STR_Create("mmore.results.mydomout");
  args->is_mytimeout = false;
  args->mytime_fileout = STR_Create("mmore.results.mytimeout");
  args->is_mythreshout = false;
  args->mythresh_fileout = STR_Create("mmore.results.mythreshout");
  // args->is_customout = false;
  // args->customout_fileout = STR_Create("results.customout");

  /* --- RANGE OPTIONS --- */
  args->t_range = (RANGE){-1, -1};
  args->q_range = (RANGE){-1, -1};
  // args->use_range = false;
  args->list_range = (RANGE){-1, -1};

  /* --- MMORE / FB-PRUNER --- */
  args->alpha = 12.0f;
  args->beta = 16.0f;
  args->gamma = 5;
  args->hard_limit = -12.0f;
  args->mmore_evalue = 2e2f;
  args->mmore_pvalue = 1e-3f;

  /* --- MMSEQS --- */
  args->mmseqs_maxhits = 1000; /* mmseqs default: 300 */
  args->mmseqs_altalis = 0;
  args->mmseqs_kmer = 7;
  args->mmseqs_kscore = 80;       /* mmseqs default: 95 */
  args->mmseqs_ungapped_vit = 0;  /* mmseqs default: 15 */
  args->mmseqs_sensitivity = 7.5; /* mmseqs default: 4 (overrides kscore) */
  args->mmseqs_p2s_pvalue = 1e-2f;
  args->mmseqs_s2s_pvalue = 1e7;
  args->mmseqs_pvalue = 1e-2f; /* mmseqs default: 1E-3 */
  args->mmseqs_evalue = 1e-2f;

  /* --- SCORE FILTERS (P_VALUES) --- */
  args->is_run_filter = false;
  args->is_run_viterbi_filter = true;
  args->is_run_cloud_filter = true;
  args->is_run_boundfwd_filter = true;
  args->is_run_fwdback_filter = true;
  args->is_run_report_filter = true;
  args->threshold_pre = 80;      /* mmseqs default: 95 */
  args->threshold_ungapped = 15; /* mmseqs default: 15 */
  args->threshold_p2s = 1e-2f;   /* mmseqs default: 1E-3 */
  args->threshold_s2s = INF;
  args->threshold_vit = 1e-2f;
  args->threshold_cloud = 1e-3f;
  args->threshold_boundfwd = 1e-4f;
  args->threshold_fwd = 1e-4f;
  args->threshold_report_eval = 10;
}

void ARGS_Dump(ARGS* args,
               FILE* fp) {
  int pad = 20;
  bool align = -1; /* -1 for right alignment, 1 for left alignment */

  fprintf(fp, "# === MMORE OPTIONS =======================\n");
  /* --- PIPELINE --- */
  fprintf(fp, "# === PIPELINE ===\n");
  fprintf(fp, "# %*s:\t%s [%d]\n", align * pad, "PIPELINE", PIPELINES[args->pipeline_mode].name, args->pipeline_mode);
  fprintf(fp, "# %*s:\t%s [%d]\n", align * pad, "VERBOSITY_MODE", VERBOSITY_NAMES[args->verbose_level], args->verbose_level);
  fprintf(fp, "# === PROGRAMS ===\n");
  fprintf(fp, "# %*s:\t%s\n", align * pad, "MMORESEQS_PROGRAM", args->mmoreseqs_program);
  fprintf(fp, "# %*s:\t%s\n", align * pad, "MMSEQS_PROGRAM", args->mmseqs_program);
  fprintf(fp, "# %*s:\t%s\n", align * pad, "HMMER_PROGRAM", args->hmmer_program);
  fprintf(fp, "# === SCRIPTS ===\n");
  fprintf(fp, "# %*s:\t%s\n", align * pad, "PREP_SCRIPT", PREP_SCRIPT);
  fprintf(fp, "# %*s:\t%s\n", align * pad, "PREPSEARCH_SCRIPT", PREPSEARCH_SCRIPT);
  fprintf(fp, "# %*s:\t%s\n", align * pad, "SEARCH_SCRIPT", SEARCH_SCRIPT);
  fprintf(fp, "# %*s:\t%s\n", align * pad, "EASYSEARCH_SCRIPT", EASYSEARCH_SCRIPT);
  fprintf(fp, "# === PIPELINE OPTIONS ===\n");
  fprintf(fp, "# %*s:\t[%d] : [%d] [%d]\n",
          align * pad, "RUN_MMSEQS", args->is_run_mmseqs, args->is_run_mmseqs_pref, args->is_run_mmseqs_align);
  fprintf(fp, "# %*s:\t[%d]\n", align * pad, "RUN_MMORE", args->is_run_mmore);
  /* --- INPUT --- */
  fprintf(fp, "# === INPUT ===\n");
  fprintf(fp, "# %*s:\t%s [%s]\n", align * pad, "TARGET", args->t_filein, FILETYPE_NAME_Get(args->t_filetype));
  fprintf(fp, "# %*s:\t%s [%s]\n", align * pad, "QUERY", args->q_filein, FILETYPE_NAME_Get(args->q_filetype));
  fprintf(fp, "# %*s:\t%s [%s]\n", align * pad, "TARGET_PREP", args->target_prep, FILETYPE_NAME_Get(args->target_prep_type));
  fprintf(fp, "# %*s:\t%s [%s]\n", align * pad, "QUERY_PREP", args->query_prep, FILETYPE_NAME_Get(args->query_prep_type));
  fprintf(fp, "# %*s:\t%s [%s]\n", align * pad, "TARGET_MMORE", args->t_mmore_filein, FILETYPE_NAME_Get(args->t_mmore_filetype));
  fprintf(fp, "# %*s:\t%s [%s]\n", align * pad, "QUERY_MMORE", args->q_mmore_filein, FILETYPE_NAME_Get(args->q_mmore_filetype));
  fprintf(fp, "# %*s:\t%s [%s]\n", align * pad, "TARGET_MMSEQS_P", args->t_mmseqs_p_filein, FILETYPE_NAME_Get(args->t_mmseqs_p_filetype));
  fprintf(fp, "# %*s:\t%s [%s]\n", align * pad, "TARGET_MMSEQS_S", args->t_mmseqs_s_filein, FILETYPE_NAME_Get(args->t_mmseqs_s_filetype));
  fprintf(fp, "# %*s:\t%s [%s]\n", align * pad, "QUERY_MMSEQS", args->q_mmseqs_filein, FILETYPE_NAME_Get(args->t_mmseqs_s_filetype));
  fprintf(fp, "# %*s:\t%s\n", align * pad, "T_INDEX_PATH", args->t_index_filein);
  fprintf(fp, "# %*s:\t%s\n", align * pad, "Q_INDEX_PATH", args->q_index_filein);
  fprintf(fp, "# %*s:\t%s\n", align * pad, "MMSEQS_M8", args->mmseqs_m8_filein);
  fprintf(fp, "# %*s:\t%s\n", align * pad, "TMP_FOLDER", args->tmp_folderpath);
  fprintf(fp, "# %*s:\t%s\n", align * pad, "PREP_FOLDER", args->prep_folderpath);
  fprintf(fp, "# \n");
  /* --- MMSEQS --- */
  fprintf(fp, "# === MMSEQS ===\n");
  fprintf(fp, "# %*s:\t%d\n", align * pad, "MMSEQS_KMER", args->mmseqs_kmer);
  fprintf(fp, "# %*s:\t%d [%d]\n", align * pad, "MMSEQS_KSCORE", args->mmseqs_kscore, !args->is_run_mmseqs_sens);
  fprintf(fp, "# %*s:\t%.3f [%d]\n", align * pad, "MMSEQS_SENS", args->mmseqs_sensitivity, args->is_run_mmseqs_sens);
  fprintf(fp, "# %*s:\t%d [%d]\n", align * pad, "MMSEQS_UNGAPPED", args->mmseqs_ungapped_vit, args->is_run_mmseqs_ungapped);
  fprintf(fp, "# %*s:\t%d\n", align * pad, "MMSEQS_ALTALIS", args->mmseqs_altalis);
  fprintf(fp, "# %*s:\t%.2e\n", align * pad, "MMSEQS_P2S_PVAL", args->threshold_p2s);
  fprintf(fp, "# %*s:\t%.2e\n", align * pad, "MMSEQS_S2S_PVAL", args->threshold_s2s);
  fprintf(fp, "# %*s:\t[%d]\n", align * pad, "MMSEQS_VITALN", args->is_run_mmseqsaln);
  /* --- MMORE / FB-PRUNER --- */
  fprintf(fp, "# === MMORE ===\n");
  fprintf(fp, "# %*s:\t%.2f\n", align * pad, "MMORE_ALPHA", args->alpha);
  fprintf(fp, "# %*s:\t%.2f\n", align * pad, "MMORE_BETA", args->beta);
  fprintf(fp, "# %*s:\t%d\n", align * pad, "MMORE_GAMMA", args->gamma);
  fprintf(fp, "# %*s:\t%.2f\n", align * pad, "MMORE_HARD_LIMIT", args->hard_limit);
  fprintf(fp, "# %*s:\t%.2e [%d]\n", align * pad, "MMORE_VITERBI_PVAL", args->threshold_vit, args->is_run_viterbi_filter);
  fprintf(fp, "# %*s:\t%.2e [%d]\n", align * pad, "MMORE_CLOUD_PVAL", args->threshold_cloud, args->is_run_cloud_filter);
  fprintf(fp, "# %*s:\t%.2e [%d]\n", align * pad, "MMORE_BOUNDFWD_PVAL", args->threshold_boundfwd, args->is_run_boundfwd_filter);
  fprintf(fp, "# %*s:\t%.2e [%d]\n", align * pad, "MMORE_REPORT_EVAL", args->threshold_report_eval, args->is_run_report_filter);

  fprintf(fp, "# %*s:\t(%d,%d)\n", align * pad, "MMORE_RANGE", args->list_range.beg, args->list_range.end);
  fprintf(fp, "# %*s:\t[%d]\n", align * pad, "MMORE_FULL", args->is_run_full);
  fprintf(fp, "# %*s:\t[%d]\n", align * pad, "MMORE_VIT_MMORE", args->is_run_vit_mmore);
  fprintf(fp, "# %*s:\t[%d]\n", align * pad, "MMORE_DOMAINS", args->is_run_domains);
  fprintf(fp, "# %*s:\t[%d] [%d]\n", align * pad, "MMORE_VITALN", args->is_run_vitaln, args->is_run_vit);
  fprintf(fp, "# %*s:\t[%d] [%d]\n", align * pad, "MMORE_POSTALN", args->is_run_postaln, args->is_run_optacc);
  fprintf(fp, "# \n");
  /* --- OUTPUT --- */
  fprintf(fp, "# === OUTPUT ===\n");
  fprintf(fp, "# %*s:\t%s [%d]\n", align * pad, "OUTPUT_FILEPATH", args->stdout_fileout, args->is_redirect_stdout);
  fprintf(fp, "# %*s:\t%s [%d]\n", align * pad, "TBLOUT_FILEPATH", args->hmmerout_fileout, args->is_hmmerout);
  fprintf(fp, "# %*s:\t%s [%d]\n", align * pad, "M8OUT_FILEPATH", args->m8out_fileout, args->is_m8out);
  fprintf(fp, "# %*s:\t%s [%d]\n", align * pad, "MYOUT_FILEPATH", args->myout_fileout, args->is_myout);
  fprintf(fp, "# %*s:\t%s [%d]\n", align * pad, "MYDOMOUT_FILEPATH", args->mydom_fileout, args->is_mydom);
  fprintf(fp, "# %*s:\t%s [%d]\n", align * pad, "MYTHRESHOUT_FILEPATH", args->mythresh_fileout, args->is_mythreshout);
  fprintf(fp, "# %*s:\t%s [%d]\n", align * pad, "MYTIMEOUT_FILEPATH", args->mytime_fileout, args->is_mytimeout);

  fprintf(fp, "# ==============================================\n\n");
}

FILETYPE ARGS_FindFiletype(STR filename) {
  for (int i = 0; i < NUM_FILETYPE_EXTS; i++) {
    char* ext_name = FILETYPE_EXTS[i].s;
    int ext_type = FILETYPE_EXTS[i].i;
    if (STRING_EndsWith(filename, ext_name, strlen(ext_name)) == 0) {
      return ext_type;
    }
  }

  fprintf(stderr, "WARNING: '%s' filetype could not be found. Set to TYPE_NULL.\n", filename);

  return FILE_NULL;
}

void ARGS_Help_Info() {
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
}

void ARGS_Command_Help_Info(ARGS* args) {
  fprintf(stdout, "USAGE: %s\n", PIPELINES_ARG_HELP[args->pipeline_mode]);
}

void ARGS_Version_Info() {
  fprintf(stdout, "MMORESEQS VERSION: v%s\n", BUILD_VERSION);
  fprintf(stdout, "BUILD_TYPE: %s\n", BUILD_TYPE);
  fprintf(stdout, "BUILD_HASH: %s\n", BUILD_HASH);

  ERRORCHECK_exit(EXIT_SUCCESS);
}

STATUS_FLAG ARGS_SetOptions(ARGS* args,
                            ARG_OPTS* arg_opts) {
  /* add all commandline options */
  PTR arg_locs[] = {&args->t_index_filein, &args->q_index_filein};
  INT arg_dtypes[] = {DATATYPE_STRING, DATATYPE_STRING};

  {
    ARG_OPTS_AddOption(arg_opts, "VERBOSE",
                       "Level of Output.",
                       "Level of Output: [0] Minimal Output [1] Errors [2] Errors+Warnings [3] Maximal Output. ( Default: [1] )",
                       "--verbose", "-v",
                       1,
                       (PTR[]){&args->verbose_level},
                       (INT[]){DATATYPE_INT});
  }
  {
    ARG_OPTS_AddOption(arg_opts, "NUM_THREADS",
                       "Number of Threads.",
                       "Level of Output: [0] Minimal Output [1] Errors [2] Errors+Warnings [3] Maximal Output. ( Default: [1] )",
                       "--num-threads", NULL,
                       1,
                       (PTR[]){&args->num_threads},
                       (INT[]){DATATYPE_INT});
  }
  { ARG_OPTS_AddOption(arg_opts, "ENFORCE WARNINGS",
                       "Should program terminate on warnings?",
                       "Should program terminate on warnings? [0] Bypass warnings [1] Terminate on warnings.  ( Default: [1] )",
                       "--enforce-warnings", "-w",
                       1,
                       (PTR[]){&args->enforce_warnings},
                       (INT[]){DATATYPE_BOOL}); }
  { ARG_OPTS_AddOption(arg_opts, "INDEX FILE",
                       "Index file location.",
                       "Index file location [query_index, target_index].  ( Default: NULL )",
                       "--index", NULL,
                       2,
                       (PTR[]){&args->t_index_filein, &args->q_index_filein},
                       (INT[]){DATATYPE_STRING, DATATYPE_STRING}); }
  { ARG_OPTS_AddOption(arg_opts, "USE LOCAL TOOLS",
                       "Should local or system tools be used?",
                       "Should local or system tools be used? [0] Use system tools [1] Use local tools.  ( Default: [1] )",
                       "--use-local-tools", NULL,
                       1,
                       (PTR[]){&args->is_use_local_tools},
                       (INT[]){DATATYPE_BOOL}); }
}
