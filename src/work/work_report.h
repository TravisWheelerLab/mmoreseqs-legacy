/*******************************************************************************
 *  - FILE:      work_report.h
 *  - DESC:    Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various
 *functions. Subroutines for outputting reports.
 *******************************************************************************/

#ifndef _WORK_REPORT
#define _WORK_REPORT

/*! FUNCTION:  	WORK_open()
 *  SYNOPSIS:  	Open all valid files in <worker>.
 */
void WORK_open(WORKER* worker);

/*! FUNCTION:  	WORK_close()
 *  SYNOPSIS:  	Close all open files in <worker>.
 */
void WORK_close(WORKER* worker);

/*! FUNCTION:  	WORK_report_header()
 *  SYNOPSIS:  	Write all report headers to all open files in <worker>.
 */
void WORK_report_header(WORKER* worker);

/*! FUNCTION:  	WORK_report_result_current()
 *  SYNOPSIS:  	Write current result entry to all open files in <worker>.
 */
void WORK_report_result_current(WORKER* worker);

/* print header for results file (default) */
/*! FUNCTION:  	WORK_report_footer()
 *  SYNOPSIS:  	Write all report footers to all open files in <worker>.
 */
void WORK_report_footer(WORKER* worker);

#endif /* _WORK_REPORT */
