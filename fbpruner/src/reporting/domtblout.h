/*******************************************************************************
 *  FILE:      domtblout.h
 *  PURPOSE:   Reporting Subroutines for generating output.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _DOMTBLOUT_H
#define _DOMTBLOUT_H

/*!   FUNCTION:   REPORT_tblout_header()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
STATUS_FLAG 
REPORT_domtblout_header(   WORKER*  worker,
                           FILE*    fp );

/*!   FUNCTION:   REPORT_alignment()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
STATUS_FLAG 
REPORT_domtblout_entry(    WORKER*  worker,
                           RESULT*  result,
                           FILE*    fp );

/*!   FUNCTION:   REPORT_tblout_header()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)              
 */
STATUS_FLAG 
REPORT_domtblout_footer(   WORKER*  worker,
                           FILE*    fp );

#endif /* _DOMTBLOUT_H */