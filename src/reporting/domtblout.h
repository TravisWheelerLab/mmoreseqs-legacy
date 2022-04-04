/*******************************************************************************
 *  - FILE:  domtblout.h
 *  - DESC:  Reporting Subroutines for generating output.
 *******************************************************************************/

#ifndef _DOMTBLOUT_H
#define _DOMTBLOUT_H

/*!   FUNCTION:   REPORT_domtblout_header()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER domtblout, see example)
 */
STATUS_FLAG
REPORT_domtblout_header(WORKER* worker, FILE* fp);

/*!   FUNCTION:   REPORT_domtblout_entry()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER domtblout, see example)
 */
STATUS_FLAG
REPORT_domtblout_entry(WORKER* worker, RESULT* result, FILE* fp);

/*!   FUNCTION:   REPORT_domtblout_footer()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER domtblout, see example)
 */
STATUS_FLAG
REPORT_domtblout_footer(WORKER* worker, FILE* fp);

#endif /* _DOMTBLOUT_H */
