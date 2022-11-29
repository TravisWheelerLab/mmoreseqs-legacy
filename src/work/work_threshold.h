/*******************************************************************************
 *  - FILE:  work_threshold.h
 *  - DESC:  Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions. Computes and tests threshold cutoffs for inner loop.
 *******************************************************************************/

#ifndef _WORK_THRESHOLD
#define _WORK_THRESHOLD

/*! FUNCTION:  	WORK_thresholds_pval_to_eval()
 *  SYNOPSIS:  	Converts threshold scores from P-values to E-values.
 */
void WORK_thresholds_pval_to_eval(WORKER* worker);

/*! FUNCTION:  	WORK_viterbi_natsc_to_eval()
 *  SYNOPSIS:  	Converts Viterbi natscore to e-value.
 */
void WORK_viterbi_mmore_natsc_to_eval(WORKER* worker);

/*! FUNCTION:  	WORK_viterbi_test_threshold()
 *  SYNOPSIS:  	Tests viterbi score against threshold.
 *  RETURN:       If passed threshold.
 */
bool WORK_viterbi_test_threshold(WORKER* worker);

/*! FUNCTION:  	WORK_cloud_natsc_to_eval()
 *  SYNOPSIS:  	Converts Cloud natscore to e-value.
 */
void WORK_cloud_natsc_to_eval(WORKER* worker);

/*! FUNCTION:  	WORK_cloud_test_threshold()
 *  SYNOPSIS:  	Tests cloud eval against threshold.
 *  RETURN:       If passed threshold.
 */
bool WORK_cloud_test_threshold(WORKER* worker);

/*! FUNCTION:  	WORK_bound_fwdback_natsc_to_eval()
 *  SYNOPSIS:  	Converts Forward natscore to e-value.
 */
void WORK_bound_fwdback_natsc_to_eval(WORKER* worker);

/*! FUNCTION:  	WORK_bound_fwdback_test_threshold()
 *  SYNOPSIS:  	Tests cloud eval against threshold.
 *  RETURN:       If passed threshold.
 */
bool WORK_bound_fwdback_test_threshold(WORKER* worker);

/*! FUNCTION:  	WORK_report_test_threshold()
 *  SYNOPSIS:  	Tests report eval against threshold (report is the
 * bias-corrected, domain-corrected score). RETURN:       <true/false> if passed
 * threshold.
 */
bool WORK_report_test_threshold(WORKER* worker);

#endif /* _WORK_THRESHOLD */
