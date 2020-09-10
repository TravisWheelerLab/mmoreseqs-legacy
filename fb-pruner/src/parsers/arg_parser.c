/*******************************************************************************
 *     FILE:   arg_parser.h
 *  PURPOSE:   Parses command line arguments. 
 *
 *  AUTHOR:    Dave Rich
 *     BUG:    
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
#include "structs.h"
#include "utilities.h"
#include "objects.h"

/* header */
#include "parsers.h"

/* Parses Arguments from the command line */
void   ARGS_Parse( ARGS*   args,
                   int     argc, 
                   char*   argv[] )
{
   int   num_main_args  = 2; 
   char* flag           = NULL;
   int   req_args       = 0;

   ARGS_Set_Defaults(args);

   /* if no <command> argument given, run test case if in debug mode */
   if (argc < 2) {
      printf("Usage: fb-pruner <command> <target_hmm_file> <query_fasta_file>\n");
      #if DEBUG 
         printf("Using DEFAULT arguments...\n\n");
         return;
      #else
         exit(EXIT_FAILURE);
      #endif
   }

   /* check for help flag */
   if ( strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0 ) {
      ARGS_Help_Info();
   }

   /* first argument is command pipeline (id_number or name) */
   bool found_pipeline = false;
   if ( isdigit(argv[1][0]) ) {
      args->pipeline_mode = atoi(argv[1]);
   } else {
      for (int i = 0; i < NUM_PIPELINE_MODES; i++) {
         if ( strcmp(argv[1], PIPELINE_NAMES[i]) == 0 ) {
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
         fprintf(stderr, "%s, ", PIPELINE_NAMES[i]);
      fprintf(stderr, "]\n");
      exit(EXIT_FAILURE);
   }

   /* set number of main arguments based on given pipeline */
   num_main_args = PIPELINE_NUM_ARGS[args->pipeline_mode];
   /* check proper number of main args remain */
   if ( argc < 2 + num_main_args ) {
      fprintf(stderr, "ERROR: Improper number of main args. [required: %d]\n", num_main_args);
      #if DEBUG 
         printf("Using DEFAULT arguments...\n\n");
         return;
      #else 
         exit(EXIT_FAILURE);
      #endif
   }

   /* number of args based on pipeline */
   if ( num_main_args == 2 )
   {
      /* second arg is query */
      ERROR_free( args->t_filepath );
      args->t_filepath = strdup(argv[2]);
      /* third arg is target */
      ERROR_free( args->q_filepath );
      args->q_filepath = strdup(argv[3]);
   }
   else if ( num_main_args == 1 )
   {
      /* second arg is query */
      args->t_filepath = strdup(argv[2]);
   }

   /* parse flags and options */
   for (int i = 2 + num_main_args; i < argc; ++i)
   {
      /* if long flag */
      if ( strncmp(argv[i], "--", 2) == 0 ) 
      {
         if ( strcmp(argv[i], (flag = "--help") ) == 0 ) {
            ARGS_Help_Info();
         }
         else if ( strcmp(argv[i], (flag = "--alpha") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->alpha = atof(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               exit(EXIT_FAILURE);
            }
         }
         else if ( strcmp(argv[i], (flag = "--beta") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->beta = atof(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               exit(EXIT_FAILURE);
            }
         }
         else if ( strcmp(argv[i], (flag = "--gamma") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->gamma = atoi(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               exit(EXIT_FAILURE);
            }
         }
         else if ( strcmp(argv[i], (flag = "--verbose") ) == 0 ) {
            req_args = 1;
            if (i+req_args < argc) {
               i++;
               args->verbose_level = atoi(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               exit(EXIT_FAILURE);
            }
         }
         else if ( strcmp(argv[i], (flag = "--eval") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               args->threshold_fbpruner = atof(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument=.\n", flag, req_args);
               exit(EXIT_FAILURE);
            }
         }
         else if ( strcmp(argv[i], (flag = "--mmseqs-tmp") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               ERROR_free(args->mmseqs_tmp_filepath);
               args->mmseqs_tmp_filepath = strdup(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument=.\n", flag, req_args);
               exit(EXIT_FAILURE);
            }
         }
         else if ( strcmp(argv[i], (flag = "--mmseqs-m8") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               ERROR_free(args->mmseqs_res_filepath);
               args->mmseqs_res_filepath = strdup(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               exit(EXIT_FAILURE);
            }
         }
         else if ( strcmp(argv[i], (flag = "--mmseqs-m8+") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               ERROR_free(args->mmseqs_plus_filepath);
               args->mmseqs_plus_filepath = strdup(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               exit(EXIT_FAILURE);
            }
         }
         else if ( strcmp(argv[i], (flag = "--mmseqs-range") ) == 0 ) {
            req_args = 2;
            if (i+req_args < argc) {
               i++;
               args->mmseqs_range.beg = atoi(argv[i]);
               i++;
               args->mmseqs_range.end = atoi(argv[i]);
            } else {
               fprintf(stderr, "ERROR: '%s' flag requires (%d) argument.\n", flag, req_args);
               exit(EXIT_FAILURE);
            }
         }
         else if ( strcmp(argv[i], (flag = "--mmseqs-lookup") ) == 0 ) {
            req_args = 2;
            if (i+req_args <= argc) {
               i++;
               ERROR_free(args->t_lookup_filepath);
               args->t_lookup_filepath = strdup(argv[i]);
               i++;
               ERROR_free(args->q_lookup_filepath);
               args->q_lookup_filepath = strdup(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               exit(EXIT_FAILURE);
            }
         }
         else if ( strcmp(argv[i], (flag = "--index") ) == 0 ) {
            req_args = 2;
            if (i+req_args <= argc) {
               i++;
               ERROR_free(args->t_indexpath);
               args->t_indexpath = strdup(argv[i]);
               i++;
               ERROR_free(args->q_indexpath);
               args->q_indexpath = strdup(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument=.\n", flag, req_args);
               exit(EXIT_FAILURE);
            }
         }
         else if ( strcmp(argv[i], (flag = "--tblout") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               ERROR_free(args->tblout_filepath);
               args->tblout_filepath = strdup(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               exit(EXIT_FAILURE);
            }
         }
         else if ( strcmp(argv[i], (flag = "--m8out") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               ERROR_free(args->m8out_filepath);
               args->m8out_filepath = strdup(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               exit(EXIT_FAILURE);
            }
         }
         else if ( strcmp(argv[i], (flag = "--myout") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               ERROR_free(args->myout_filepath);
               args->myout_filepath = strdup(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               exit(EXIT_FAILURE);
            }
         }
         else if ( strcmp(argv[i], (flag = "--output") ) == 0 ) {
            req_args = 1;
            if (i+req_args <= argc) {
               i++;
               args->output_filepath = strdup(argv[i]);
            } else {
               fprintf(stderr, "ERROR: %s flag requires (%d) argument.\n", flag, req_args);
               exit(EXIT_FAILURE);
            }
         }
         else {
            fprintf(stderr, "ERROR: '%s' is not a recognized flag.\n", argv[i]);
            exit(EXIT_FAILURE);
         }
      }
      else
      {
         fprintf(stderr, "ERROR: '%s' is not associated with an argument or flag\n", argv[i]);
         exit(EXIT_FAILURE);
      }
   }

   args->t_filetype  = ARGS_Find_FileType( args->t_filepath );
   args->q_filetype  = ARGS_Find_FileType( args->q_filepath );
}

/* SET DEFAULT ARGUMENTS (generic) */
void  ARGS_Set_Defaults( ARGS* args )
{
   args->qt_search_space         = SELECT_ALL_V_ALL;

   args->t_filepath              = strdup("test-input/test1_2.hmm");
   args->q_filepath              = strdup("test-input/test1_1.fa");

   args->t_indexpath             = NULL;
   args->q_indexpath             = NULL;

   args->t_filetype              = FILE_HMM;
   args->q_filetype              = FILE_FASTA;

   /* TODO: tblout testing */
   args->is_redirect_stdout      = false;
   args->output_filepath         = strdup("results.stdout");

   args->is_tblout               = false;
   args->tblout_filepath         = strdup("results.tblout");

   args->is_m8out                = true;
   args->m8out_filepath          = strdup("results.m8");

   args->is_myout                = false;
   args->myout_filepath          = strdup("results.myout");

   args->mmseqs_res_filepath     = NULL;
   args->mmseqs_plus_filepath    = NULL;
   args->mmseqs_tmp_filepath     = NULL;

   args->tmp_folderpath          = NULL;

   args->tmp_remove              = false;
   args->dbg_folderpath          = strdup("test-output/");

   args->alpha                   = 12.0f;
   args->beta                    = 20.0f;
   args->gamma                   = 5;

   args->pipeline_mode           = PIPELINE_TEST;
   args->verbose_level           = VERBOSE_LOW;
   args->search_mode             = MODE_UNILOCAL;

   args->filter_on               = false;
   args->threshold_vit           = 1e-3f;
   args->threshold_fwd           = 1e-5f;
   args->threshold_fbpruner      = 1e-5f;

   /* these will default to entire file unless filled with positive ints */
   args->t_range                 = (RANGE) { -1, -1 };    
   args->q_range                 = (RANGE) { -1, -1 };
   args->mmseqs_range            = (RANGE) { -1, -1 };
}


/* sends ARGS data to FILE POINTER */
void ARGS_Dump( ARGS*    args,
                FILE*    fp )
{
   int      pad               = 20;
   bool     align             = 1;     /* -1 for right alignment, 1 for left alignment */

   fprintf( fp, "=== ARGS =====================\n");
   /* parameters */
   fprintf( fp, "%*s:\t%s == %d\n",    align * pad,  "PIPELINE",        PIPELINE_NAMES[args->pipeline_mode],   args->pipeline_mode );
   fprintf( fp, "%*s:\t%s == %d\n",    align * pad,  "VERBOSITY_MODE",  VERBOSITY_NAMES[args->verbose_level],  args->verbose_level );
   fprintf( fp, "%*s:\t%s\n",          align * pad,  "SEARCH_MODE",     MODE_NAMES[args->search_mode] );
   fprintf( fp, "%*s:\t%.3f\n",        align * pad,  "ALPHA",           args->alpha );
   fprintf( fp, "%*s:\t%.3f\n",        align * pad,  "BETA",            args->beta );
   fprintf( fp, "%*s:\t%d\n",          align * pad,  "GAMMA",           args->gamma );
   fprintf( fp, "\n" );
   /* inputs */
   fprintf( fp, "%*s:\t%s\n",          align * pad,  "TARGET_FILEPATH", args->t_filepath );
   fprintf( fp, "%*s:\t%s\n",          align * pad,  "TARGET_FILETYPE", FILE_TYPE_NAMES[args->t_filetype] );
   fprintf( fp, "%*s:\t%s\n",          align * pad,  "QUERY_FILEPATH",  args->q_filepath );
   fprintf( fp, "%*s:\t%s\n",          align * pad,  "QUERY_FILETYPE",  FILE_TYPE_NAMES[args->q_filetype] );
   fprintf( fp, "\n" );
   /* index input */
   fprintf( fp, "%*s:\t%s\n",          align * pad,  "T_INDEX_PATH",    args->t_indexpath );
   fprintf( fp, "%*s:\t%s\n",          align * pad,  "Q_INDEX_PATH",    args->q_indexpath );
   fprintf( fp, "\n" );
   /* mmseqs input */
   fprintf( fp, "%*s:\t%s\n",          align * pad,  "MMSEQS_RESULTS",  args->mmseqs_res_filepath );
   fprintf( fp, "%*s:\t%s\n",          align * pad,  "MMSEQS_TEMP",     args->mmseqs_tmp_filepath );  
   fprintf( fp, "\n" );
   /* output */
   fprintf( fp, "%*s:\t%s\n",          align * pad,  "OUTPUT_FILEPATH", args->output_filepath );
   fprintf( fp, "=============================\n\n");
}

/* examines target and query, and finds the type of the files */
int ARGS_Find_FileType( char* _filename_ )
{
   for (int i = 0; i < NUM_FILE_EXTS; i++) {
      char* ext = FILE_TYPE_EXTS[i];
      if ( STRING_EndsWith( _filename_, ext, strlen(ext) ) == 0 ) {
         return FILE_TYPE_MAP[i];
      }
   }

   fprintf(stderr, "ERROR: '%s' is not an acceptable file type.\n", _filename_);
   exit(EXIT_FAILURE);
   return -1;
}

/* output help info */
void ARGS_Help_Info()
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
   exit(EXIT_SUCCESS);
}

/* output version info */
void ARGS_Version_Info()
{
   
}