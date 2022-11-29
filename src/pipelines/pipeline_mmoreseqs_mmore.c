/*******************************************************************************
 *  - FILE:   pipeline_mmseqs.c
 *  - DESC:    MMORE-SEQS core pipeline.
 *          Post-MMSEQS phase, it takes in the results of that pipeline.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <stdbool.h>

/* local imports */
#include "../objects/structs.h"
#include "../work/_work.h"

/* header */
#include "_pipelines.h"


/* private functions */
STATUS_FLAG mmore_main_SetDefault_Tasks(TASKS* tasks);


STATUS_FLAG mmoreseqs_mmore_pipeline(WORKER* worker) {
  ARGS* args = worker->args;

  printf("=== MMORESEQS: MMORE SEARCH PIPELINE ===\n");

  WORK_init(worker);

  WORK_open(worker);

  mmore_main_SetDefault_Tasks(worker->tasks);

  WORK_load_indexes(worker);

  // update thresholds to e-values to reflect database size
  WORK_thresholds_pval_to_eval(worker);

  WORK_load_mmseqs_file(worker);

  WORK_report_header(worker);

  // get bounds of mmseqs results to be searched
  int i_cnt = 0;
  int i_beg = args->list_range.beg;
  int i_end = args->list_range.end;
  int i_rng = i_end - i_beg;

  args->verbose_level = 1000;

  printf("# Beginning search through mmseqs-m8 list on range (%d,%d)...\n", i_beg, i_end);

  // threshold tests
  bool passed[4];

  /* === ITERATE OVER EACH RESULT === */
  /* Look through each input result (i = index in full list, i_cnt = index relative to search range) */
  for (int i = i_beg; i < i_end; i++, i_cnt++) {
    printf(
      "\n# (%d/%d): Running cloud search for result (%d of %d)...\n",
      i_cnt, i_rng, i + 1, i_end
    );

    passed[0] = false;
    passed[1] = false;
    passed[2] = false;
    passed[3] = false;

    /* prep for current iteration */
    WORK_preiter(worker);

    /* get next mmseqs entry */
    WORK_load_mmseqs_by_id(worker, i_cnt);

    /* evaluate mmseqs viterbi scoring filter */
    passed[0] = WORK_viterbi_test_threshold(worker);

    /* check if mmseqs viterbi passes threshold */
    if (passed[0] == true) {
      printf(":: VITERBI PASSED ::\n");
      /* load target hmm profile from file */
      WORK_load_target(worker);

      /* load query sequence from file */
      WORK_load_query(worker);

      /* clear old data and update data structs for problem size */
      WORK_reuse(worker);

      /* get viterbi alignment bounds from mmseqs entry */
      WORK_load_mmseqs_alignment(worker);

      /* run cloud search */
      WORK_cloud_search_linear(worker);

      /* evaluate cloud search scoring filter */
      WORK_cloud_natsc_to_eval(worker);

      passed[1] = WORK_cloud_test_threshold(worker);
    }

    /* extra work */
    if (args->is_run_vit_mmore == true) {
      worker->tasks->lin_vit = true;
      WORK_viterbi_mmore(worker);
      /* evaluate viterbi */
      WORK_viterbi_mmore_natsc_to_eval(worker);
    }

    /* check if cloud search composite score passes threshold */
    if (passed[0] == true && passed[1] == true) {
      fprintf_vall(stdout, ":: CLOUD PASSED ::\n");
      /* merge and reorient cloud */
      WORK_cloud_merge_and_reorient(worker);
      /* run bound forward */
      WORK_bound_fwdback_linear(worker);
      /* bound forward scoring filter */
      WORK_bound_fwdback_natsc_to_eval(worker);

      passed[2] = WORK_bound_fwdback_test_threshold(worker);
    }

    /* check if bound forward score passes threshold */
    if (passed[0] == true && passed[1] == true && passed[2] == true) {
      fprintf_vall(stdout, ":: FORWARD PASSED ::\n");
      /* compute posterior and bias, find domains, compute domain-specific posterior and bias */
      WORK_posterior(worker);
      /* run posterior for each found domain */
      WORK_posterior_bydom(worker);
      /* build final scores */
      WORK_construct_scores(worker);

      passed[3] = WORK_report_test_threshold(worker);
    }

    /* print thresholds which passed */
    if (worker->args->verbose_level >= VERBOSE_HIGH) {
      printf_vall("THRESHOLDS PASSED: %d => %d => %d => %d\n",
                  passed[0], passed[1], passed[2], passed[3]);
    }

    /* cleanup for current iteration */
    WORK_postiter(worker);

    /* only report if all thresholds passed */
    if ((passed[0] == true && passed[1] == true && passed[2] == true && passed[3] == true)) {
      fprintf_vall(stdout, ":: REPORT PASSED ::\n");
      /* print results */
      WORK_report_result_current(worker);
    }
  }

  /* cleanup for end of loop */
  WORK_postloop(worker);
  /* add footer to all reports */
  WORK_report_footer(worker);
  /* close all file pointers */
  WORK_close(worker);
  /* free pipeline-specific worker data structures */
  WORK_cleanup(worker);
}

/*! FUNCTION:  	mmore_searchmmore_pipeline()
 *  SYNOPSIS:  	Pipeline for tail of MMore, post-MMseqs.
 *                Runs Adaptive Banding / Cloud Search step of MMORE pipeline.
 */
STATUS_FLAG
mmore_main_SetDefault_Tasks(TASKS* tasks) {
  /* set flags for pipeline tasks */
  /* TASKS */
  {
    /* sparse algs */
    tasks->sparse = true;
    tasks->sparse_bound_fwd = true;
    tasks->sparse_bound_bck = true;
    tasks->sparse_bias_corr = true;
    /* linear algs */
    tasks->linear = true;         /* if any other linear tasks are flagged, this must be too */
    tasks->lin_fwd = false;       /* optional, but can't recover alignment */
    tasks->lin_bck = false;       /* optional, but can't recover alignment */
    tasks->lin_vit = false;       /* optional, but can't recover alignment */
    tasks->lin_trace = false;     /* optional, but can't recover alignment */
    tasks->lin_cloud_fwd = true;  /* required for sparse and linear bound fwdbck */
    tasks->lin_cloud_bck = true;  /* required for sparse and linear bound fwdbck */
    tasks->lin_bound_fwd = true;  /* can't be used to recover alignment */
    tasks->lin_bound_bck = false; /* can't be used to recover alignment */
    /* quadratic algs */
    tasks->quadratic = false;      /* if any other quadratic tasks are flagged, this must be too */
    tasks->quad_fwd = false;       /* optional */
    tasks->quad_bck = false;       /* optional */
    tasks->quad_vit = false;       /* viterbi required for cloud search */
    tasks->quad_trace = false;     /* traceback required for cloud search  */
    tasks->quad_bound_fwd = false; /* required step of cloud search */
    tasks->quad_bound_bck = false; /* optional */
    tasks->quad_bias_corr = false; /* optional: requires quadratic forward backward */
  }
}