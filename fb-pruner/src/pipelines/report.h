/*******************************************************************************
 *  FILE:      report.h
 *  PURPOSE:   Reporting Subroutines
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _REPORT_H
#define _REPORT_H

/* print header to console */
void REPORTER_header(   WORKER*  worker,
                        FILE*    fp );

/* print horizontal rule */
void REPORTER_hr( FILE* fp );

/* print summary statistics */
void REPORTER_summary_stats(  WORKER*  worker,
                              FILE*    fp );

/* print alignment */
void REPORTER_alignment(   WORKER*  worker,
                           FILE*    fp );

#endif /* _REPORT_H */