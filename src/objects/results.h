/*******************************************************************************
 *  FILE:      results.h
 *  PURPOSE:   RESULTS object
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _RESULTS_H
#define _RESULTS_H

/* constructor */
RESULTS* RESULTS_Create();
/* destructor */
void RESULTS_Destroy( RESULTS* res );

/* add result to list */
void RESULTS_PushBack( RESULTS*  res, 
                       RESULT*   r );
/* resize results */
void RESULTS_Resize( RESULTS* res,
                     int      size );

/* output results as .m8 format to file pointer */
void RESULTS_M8_Dump( RESULTS*   res,
                      FILE*      fp );

/* output results in personal format to file pointer */
void RESULTS_My_Dump( RESULTS*   res,
                      FILE*      fp );

/* output single results to file pointer */
void RESULT_M8_Dump( RESULT*  res,
                     FILE*    fp );

#endif /* _RESULTS_H */