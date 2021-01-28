/*******************************************************************************
 *  FILE:      mythreshout.h
 *  PURPOSE:   Reporting for generating mythreshout format output.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _MYTHRESHOUT_H
#define _MYTHRESHOUT_H

/*    FUNCTION:   REPORT_mythreshout_header()
 *    SYNOPSIS:   Print all alignment data for current search
 */
void REPORT_mythreshout_header(  WORKER*  worker,
                                 FILE*    fp );

/*    FUNCTION:   REPORT_mythreshout_entry()
 */
void REPORT_mythreshout_entry(   WORKER*  worker,
                                 RESULT*  result,
                                 FILE*    fp );

/*    FUNCTION:   REPORT_mythreshout_footer()
 *    SYNOPSIS:   Print footer
 *                (modeled after HMMER, see example)              
 */
void REPORT_mythreshout_footer(  WORKER*  worker,
                                 FILE*    fp );

#endif /* _MYTHRESHOUT_H */