/*******************************************************************************
 *  - FILE:      mytimeout.h
 *  - DESC:    Reporting Subroutines for generating mytimeout.
 *******************************************************************************/

#ifndef _MYTIMEOUT_H
#define _MYTIMEOUT_H

/* === MYTIMEOUT FUNCTIONS === */

/*!  FUNCTION:    REPORT_mytimeout_header()
 *   SYNOPSIS:    Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
STATUS_FLAG
REPORT_mytimeout_header(WORKER* worker, FILE* fp);

/*!  FUNCTION:    REPORT_mytimeout_entry()
 *   SYNOPSIS:    Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
STATUS_FLAG
REPORT_mytimeout_entry(WORKER* worker, RESULT* result, FILE* fp);

/*!  FUNCTION:    REPORT_mytimeout_footer()
 *   SYNOPSIS:    Print footer
 *                (modeled after HMMER, see example)
 */
STATUS_FLAG
REPORT_mytimeout_footer(WORKER* worker, FILE* fp);

/*!  FUNCTION:    REPORT_myout_entry_multiline()
 *   SYNOPSIS:    Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
STATUS_FLAG
REPORT_mytimeout_totals(WORKER* worker, RESULT* result, FILE* fp);

#endif /* _MYTIMEOUT_H */
