/*******************************************************************************
 *  FILE:      main.c
 *  PURPOSE:   Main Method, Parses Command Line Arguments, then
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

/* data structures and utility functions */
#include "objects/structs.h"
#include "utility.h"

/* file parsers */
#include "hmm_parser.h"
#include "seq_parser.h"
#include "index_parser.h"

/* objects */
#include "objects/alignment.h"
#include "objects/sequence.h"
#include "objects/hmm_profile.h"
#include "objects/edgebound.h"
#include "objects/clock.h"
#include "objects/f_index.h"
#include "objects/mystring.h"

#include "objects/matrix/matrix_2d.h"
#include "objects/matrix/matrix_3d.h"
#include "objects/vectors/vector_range_2d.h"

/* viterbi & fwdbck (quadratic) */
#include "viterbi_quad.h"
#include "traceback_quad.h"
#include "fwdback_quad.h"

/* cloud search (naive) */
#include "bounded_fwdbck_naive.h"
/* cloud search (quadratic space) */
#include "cloud_search_quad.h"
#include "merge_reorient_quad.h"
#include "bounded_fwdbck_quad.h"
/* cloud search (linear space) */
#include "cloud_search_linear.h"
#include "merge_reorient_linear.h"
#include "bounded_fwdbck_linear.h"
/* temp test */
#include "cloud_search_linear_rows.h"

/* debugging methods */
#include "testing.h"

/* pipelines */
#include "pipeline.h"

/* set debug macros */
#ifndef DEBUG
   #define DEBUG false
#endif

/* === HEADER === */
/* Parses Arguments from the command line */
ARGS* ARGS_Parse(int     argc, 
                 char*   argv[]);
/* SET DEFAULT ARGUMENTS (for testing) */
void  ARGS_Set_Defaults(ARGS* args);
/* sends ARGS data to FILE POINTER */
void ARGS_Dump(ARGS*    args,
               FILE*    fp);
/* examines target and query, and finds the type of the files */
int determine_FileType( char* _filename_ );

/* === MAIN ENTRY-POINT TO PROGRAM === */
int main ( int argc, char *argv[] )
{
   /* parse command line arguments */
   ARGS* args  = NULL;
   args = ARGS_Parse(argc, argv);

   /* output arguments */
   ARGS_Dump( args, stdout );

   /* jumps to pipeline based on -p flag */
   printf("> Running %s...\n\n", PIPELINE_NAMES[args->pipeline_mode] );
   PIPELINES[ args->pipeline_mode ]( args );

   exit(EXIT_SUCCESS);
}

/* Parses Arguments from the command line */
ARGS*  ARGS_Parse( int     argc, 
                   char*   argv[] )
{
   ARGS*       args           = NULL;
   int         num_main_args  = 0; 
   const int   max_main_args  = 2;

   args = (ARGS*) malloc( sizeof(ARGS) );
   if (args == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc ARGS.\n");
      exit(EXIT_FAILURE);
   }
   ARGS_Set_Defaults(args);

   if (argc == 1) {
      printf("Usage: ./fb-pruner <target_hmm_file> <query_fasta_file>\n");
      printf("Using DEFAULT arguments...\n\n");
   }

   for (int i = 1; i < argc; ++i)
   {
      if ( argv[i][0] == '-' )
      {
         switch (argv[i][1]) {

            /* alpha value */
            case 'a':   
               i++;
               if (i < argc) {
                  args->alpha = atof(argv[i]);
               } else {
                  fprintf(stderr, "ERROR: -a flag requires argument: <alpha>.\n");
               }
               break;

            /* beta value */
            case 'b':  
               i++;
               if (i < argc) {
                  args->beta = atoi(argv[i]);
               } else {
                  fprintf(stderr, "ERROR: -b flag requires argument: <beta>.\n");
               }
               break;

            /* append results to outfile */
            case 'o':  
               i++;
               if (i < argc) {
                  args->output_filepath = strdup(argv[i]);
                  // args->outfile = fopen(argv[i], "w+");
               } else {
                  fprintf(stderr, "ERROR: -o flag requires argument: <output_filepath>.\n");
               } 
               break;

            /* set windows (viterbi alignments) */
            case 'w': 
               i ++;
               if (i+3 < argc) {
                  args->beg.j = atoi(argv[i]);
                  i++;
                  args->beg.i = atoi(argv[i]);
                  i++;
                  args->end.j = atoi(argv[i]);
                  i++;
                  args->end.i = atoi(argv[i]);
                  i++;
               } else {
                  fprintf(stderr, "ERROR: -w flag requires 4 arguments: <beg_j, beg_i, end_j, end_i>.\n");
               }

            /* set pipeline */
            case 'p': 
               i++;
               if (i < argc) {
                  args->pipeline_mode = atoi(argv[i]);
               } else {
                  fprintf(stderr, "ERROR: -p flag requires argument: <pipeline_mode>.\n");
               }
               break;

            /* input results file */
            case 'i':
               i++;
               if (i < argc) {
                  args->hits_filepath = strdup(argv[i]);
               } else {
                  fprintf(stderr, "ERROR: -i flag requires argument: <m8_hits_file>.\n");
               }
               break;

            /* help */
            case 'h':
               fprintf(stdout, "USAGE: ./cloud_fwdbck <target_file> <query_file>\n");
               exit(EXIT_SUCCESS);
               break;

            /* "--" options */
            case '-':
               break;

            default:
               fprintf(stderr, "ERROR: -%c is not a recognized flag.\n", argv[i][1]);
               exit(EXIT_FAILURE);
         }
      }
      else
      {
         printf("%s\n", argv[i]);
         if (num_main_args == 0)
         {
            args->target_filepath = strdup(argv[i]);
         }
         else if (num_main_args == 1)
         {
            args->query_filepath = strdup(argv[i]);
         }
         else
         {
            fprintf(stderr, "ERROR: Too many main arguments.\n");
         }
         num_main_args++;
      }
   }

   args->target_filetype = determine_FileType( args->target_filepath );
   args->query_filetype  = determine_FileType( args->query_filepath );

   return args;
}

/* SET DEFAULT ARGUMENTS (for testing) */
void  ARGS_Set_Defaults( ARGS* args )
{
   args->target_filepath         = "data/test1_2.hmm";
   args->query_filepath          = "data/test1_1.fa";

   args->target_indexpath        = NULL;
   args->query_indexpath         = NULL;

   args->target_filetype         = FILE_HMM;
   args->query_filetype          = FILE_FASTA;

   args->hits_filepath           = NULL;

   args->output_filepath         = "#stdout";
   // args->output_filepath         = "stats/cloud_stats.tsv";

   args->alpha                   = 20.0f;
   args->beta                    = 5;

   args->pipeline_mode           = PIPELINE_TIME;
   args->verbosity_mode          = VERBOSE_NONE;
   args->search_mode             = MODE_UNILOCAL;

   args->viterbi_sc_threshold    = 0.0f;
   args->fwdbck_sc_threshold     = 0.0f;
   args->cloud_sc_threshold      = 0.0f;
}

/* sends ARGS data to FILE POINTER */
void ARGS_Dump( ARGS*    args,
                FILE*    fp )
{
   int      pipeline          = args->pipeline_mode;
   int      verbosity         = args->verbosity_mode;
   int      search_mode       = args->search_mode;

   float    alpha             = args->alpha;
   int      beta              = args->beta;

   char*    t_filepath        = args->target_filepath;
   char*    q_filepath        = args->query_filepath;
   int      t_filetype        = args->target_filetype;
   int      q_filetype        = args->query_filetype;

   char*    hits_filepath     = args->hits_filepath;

   char*    output_filepath   = args->output_filepath;

   int      pad               = 20;
   bool     left_aln          = 1;

   fprintf( fp, "=== ARGS =====================\n");

   fprintf( fp, "%*s:\t%s\n",    left_aln * pad,  "PIPELINE",         PIPELINE_NAMES[pipeline] );
   fprintf( fp, "%*s:\t%s\n",    left_aln * pad,  "VERBOSITY_MODE",   VERBOSITY_NAMES[verbosity] );
   fprintf( fp, "%*s:\t%s\n",    left_aln * pad,  "SEARCH_MODE",      MODE_NAMES[search_mode] );
   fprintf( fp, "%*s:\t%.3f\n",  left_aln * pad,  "ALPHA",            alpha );
   fprintf( fp, "%*s:\t%d\n",    left_aln * pad,  "BETA",             beta );
   fprintf( fp, "\n" );

   fprintf( fp, "%*s:\t%s\n",    left_aln * pad,  "TARGET_FILEPATH",  t_filepath );
   fprintf( fp, "%*s:\t%s\n",    left_aln * pad,  "TARGET_FILETYPE",  FILE_TYPE_NAMES[t_filetype] );
   fprintf( fp, "%*s:\t%s\n",    left_aln * pad,  "QUERY_FILENAME",   q_filepath );
   fprintf( fp, "%*s:\t%s\n",    left_aln * pad,  "QUERY_FILETYPE",   FILE_TYPE_NAMES[q_filetype] );
   fprintf( fp, "%*s:\t%s\n",    left_aln * pad,  "HITS_FILEPATH",    hits_filepath );
   fprintf( fp, "\n" );

   fprintf( fp, "%*s:\t%s\n",    left_aln * pad,  "OUTPUT_FILEPATH",  output_filepath );

   fprintf( fp, "=============================\n\n");
}

/* examines target and query, and finds the type of the files */
int determine_FileType( char* _filename_ )
{
   for (int i = 0; i < NUM_FILE_EXTS; i++) {
      char* ext = FILE_TYPE_EXTS[i];
      if ( STRING_EndsWith( _filename_, ext, strlen(ext) ) == 0 ) {
         return FILE_TYPE_MAP[i];
      }
   }

   fprintf(stderr, "ERROR: '%s' is not an acceptable file type.\n", _filename_);
   exit(EXIT_FAILURE);

   return 0;
}
