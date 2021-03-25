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
RESULTS* 
RESULTS_Create();

/* destructor */
RESULTS* 
RESULTS_Destroy( RESULTS* res );

/* add result to list */
void 
RESULTS_Pushback(    RESULTS*  res, 
                     RESULT*   r );
/* resize results */
void 
RESULTS_Resize(   RESULTS* 	res,
                  size_t   	size );

/* output results in personal format to file pointer */
void 
RESULTS_Dump(     RESULTS*   res,
                  FILE*      fp );

/* output single results to file pointer */
void 
RESULT_Dump(   RESULT*  res,
               FILE*    fp );

#endif /* _RESULTS_H */