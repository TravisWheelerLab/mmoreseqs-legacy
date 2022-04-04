/*******************************************************************************
 *  - FILE:      work.c
 *  - DESC:    Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various
 *functions.
 *******************************************************************************/

#ifndef _WORK_H
#define _WORK_H

/* maintenance and outside main loop */
#include "work_index.h"
#include "work_maintenance.h"
/* main loop functions */
#include "work_loader.h"
#include "work_loop.h"
// #include "work_target_hmmprof.h"
// #include "work_query_seq.h"
#include "work_report.h"
#include "work_threshold.h"
/* algorithms */
#include "work_cloud_fwdback.h"
#include "work_cloud_merge.h"
#include "work_cloud_search.h"
#include "work_fwdback.h"
#include "work_optacc.h"
#include "work_posterior.h"
#include "work_posterior_bydom.h"
#include "work_scoring.h"
#include "work_viterbi.h"
/* scripting */
#include "work_scripting.h"
/* other */
#include "work_etc.h"
#include "work_sparse_mx.h"

#endif /* _WORK_H */
