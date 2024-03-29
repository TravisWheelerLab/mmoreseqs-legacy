/*******************************************************************************
 *  - FILE:   m8_parser.c
 *  - DESC:    Parser for .m8 results file.
 *******************************************************************************/

/* imports */
#define _GNU_SOURCE
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

/* === M8 FORMAT ======================================
           qseqid means Query Seq-id
           sseqid means Subject Seq-id
           pident means Percentage of identical matches
           length means Alignment length
         mismatch means Number of mismatches
          gapopen means Number of gap openings
           qstart means Start of alignment in query
             qend means End of alignment in query
           sstart means Start of alignment in subject
             send means End of alignment in subject
           evalue means Expect value
         bitscore means Bit score
   ====================================================
 */

void RESULTS_M8_Parse(M8_RESULTS* results,
                      char* filename,
                      int start_idx,
                      int end_idx) {
  /* parser vars */
  FILE* fp = NULL;
  int line_count = -1;      /* line counter of current line in file */
  int result_count = -1;    /* tracks index of current result */
  char* line_buf = NULL;    /* pointer to start of buffered line */
  size_t line_buf_size = 0; /* length of entire <line_buf> array */
  size_t line_size = 0;     /* length of current line in <line_buf> array */

  M8_RESULT res_tmp;     /* temporary result for storing current line */
  char* line_ptr = NULL; /* moving pointer for iterating over tokens in <line_buf> */
  char* token = NULL;    /* token that tracks each word in <line_buf> */

  /* open file */
  fp = fopen(filename, "r");
  /* check for file read error */
  if (fp == NULL) {
    fprintf(stderr, "ERROR: Bad FILE POINTER for .M8 PARSER => %s\n", filename);
    ERRORCHECK_exit(EXIT_FAILURE);
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

  while ((line_size = getline(&line_buf, &line_buf_size, fp)), line_size != -1) {
    line_count++;

    /* ignore comment lines */
    if (line_buf[0] == '#' || line_size <= 1) {
      continue;
    }

    /* if we are before start or end of target results */
    result_count++;
    if (result_count < start_idx) {
      continue;
    }
    if (result_count > end_idx) {
      break;
    }

    /* remove newline from end of line */
    if (line_size > 0 && line_buf[line_size - 1] == '\n') {
      line_buf[--line_size] = '\0';
    }

    // fprintf(stdout, "[%d] %s\n", line_count, line_buf);
    res_tmp.result_id = result_count;

    /* split line on spaces, tabs, and newlines */
    line_ptr = line_buf;

    /* [1] query */
    token = strtok_r(line_ptr, " \t", &line_ptr);
    // if (res_tmp.query_name != NULL) ERROR_free(res_tmp.query_name);
    res_tmp.query_name = STR_Create(token);

    /* [2] target */
    token = strtok_r(line_ptr, " \t", &line_ptr);
    // if (res_tmp.target_name != NULL) ERROR_free(res_tmp.target_name);
    res_tmp.target_name = STR_Create(token);

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
    res_tmp.q_beg = atoi(token);

    /* [8] query end */
    token = strtok_r(line_ptr, " \t", &line_ptr);
    res_tmp.q_end = atoi(token);

    /* [9] target start */
    token = strtok_r(line_ptr, " \t", &line_ptr);
    res_tmp.t_beg = atoi(token);

    /* [10] target end */
    token = strtok_r(line_ptr, " \t", &line_ptr);
    res_tmp.t_end = atoi(token);

    /* [11] E-value */
    token = strtok_r(line_ptr, " \t", &line_ptr);
    res_tmp.eval = atof(token);

    /* [12] bit-score */
    token = strtok_r(line_ptr, " \t", &line_ptr);
    res_tmp.bitsc = atoi(token);

    /* add new result to results list */
    M8_RESULTS_Pushback(results, &res_tmp);
  }

  results->num_hits = line_count + 1;
  results->num_searches = line_count + 1;
}
