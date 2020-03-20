/*******************************************************************************
 *  FILE:      m8_parser.c
 *  PURPOSE:   Parser for .m8 results file.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

/* local imports */
#include "../structs.h"
#include "objects/results.h"

/* header */
#include "m8_parser.h"

/* parse .m8 results file and create RESULTS object */
RESULTS* RESULTS_M8_Parse( char* _filename_ )
{
   RESULTS* results  = NULL;
   results           = RESULTS_Create();

   /* parser vars */
   FILE*    fp             = NULL;    
   int      line_count     = -1;       /* line counter of current line in file */
   char*    line_buf       = NULL;     /* pointer to start of buffered line */
   size_t   line_buf_size  = 0;        /* length of entire <line_buf> array */
   size_t   line_size      = 0;        /* length of current line in <line_buf> array */

   RESULT   res_tmp;                   /* temporary result for storing current line */
   char*    line_ptr       = NULL;     /* moving pointer for iterating over tokens in <line_buf> */
   char*    token          = NULL;     /* token that tracks each word in <line_buf> */

   /* open file */
   fp = fopen(_filename_, "r");
   /* check for file read error */
   if (fp == NULL)
   {
      fprintf(stderr, "ERROR: Bad FILE POINTER for .M8 PARSER => %s\n", _filename_ );
      exit(EXIT_FAILURE);
   }

   /* initialize temporary vars */
   // res_tmp.target_name  = NULL;
   // res_tmp.query_name   = NULL;
   // res_tmp.perc_id      = 0.f;
   // res_tmp.aln_len      = 0;
   // res_tmp.mismatch     = 0;
   // res_tmp.aln_len      = 0;
   // res_tmp.gap_openings = 0.f;
   // res_tmp.query_start  = 0;
   // res_tmp.query_end    = 0;
   // res_tmp.target_start = 0;
   // res_tmp.target_end   = 0;
   // res_tmp.e_value      = 0.f;
   // res_tmp.bit_score    = 0;


   while ( ( line_size = getline ( &line_buf, &line_buf_size, fp ) ), line_size != -1 )
   {
      line_count++;

      /* ignore comment lines */
      if (line_buf[0] == '#') {
         continue;
      }
      /* remove newline from end of line */
      if (line_size > 0 && line_buf[line_size-1] == '\n') {
         line_buf[--line_size] = '\0';
      }

      // fprintf(stdout, "[%d] %s\n", line_count, line_buf);

      /* split line on spaces, tabs, and newlines */
      line_ptr = line_buf;

      /* [1] query */ 
      token = strtok_r(line_ptr, " \t", &line_ptr);
      // if (res_tmp.query_name != NULL) free(res_tmp.query_name);
      res_tmp.query_name = strdup(token);

      /* [2] target */
      token = strtok_r(line_ptr, " \t", &line_ptr);
      // if (res_tmp.target_name != NULL) free(res_tmp.target_name);
      res_tmp.target_name = strdup(token);

      /* [3] percent id */
      token = strtok_r(line_ptr, " \t", &line_ptr);
      res_tmp.perc_id = atof(token);

      /* [4] alignment length */
      token = strtok_r(line_ptr, " \t", &line_ptr);
      res_tmp.aln_len = atoi(token);

      /* [5] alignment length */
      token = strtok_r(line_ptr, " \t", &line_ptr);
      res_tmp.mismatch = atoi(token);

      /* [6] gap openings */
      token = strtok_r(line_ptr, " \t", &line_ptr);
      res_tmp.gap_openings = atoi(token);

      /* [7] query start */
      token = strtok_r(line_ptr, " \t", &line_ptr);
      res_tmp.query_start = atoi(token);

      /* [8] query end */
      token = strtok_r(line_ptr, " \t", &line_ptr);
      res_tmp.query_end = atoi(token);

      /* [9] target start */
      token = strtok_r(line_ptr, " \t", &line_ptr);
      res_tmp.target_start = atoi(token);

      /* [10] target end */
      token = strtok_r(line_ptr, " \t", &line_ptr);
      res_tmp.target_end = atoi(token);

      /* [11] E-value */
      token = strtok_r(line_ptr, " \t", &line_ptr);
      res_tmp.e_value = atof(token);

      /* [12] bit-score */
      token = strtok_r(line_ptr, " \t", &line_ptr);
      res_tmp.bit_score = atoi(token);

      /* add new result to results list */
      RESULTS_PushBack( results, &res_tmp );
   }

   return results;
}