/*******************************************************************************
 *  FILE:      m8_results.h
 *  PURPOSE:   M8_RESULTS object.
 *             Stores results in .m8 format.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _M8_RESULTS_H
#define _M8_RESULTS_H

/*! FUNCTION:  M8_RESULTS_Create()
 *  SYNOPSIS:  Create <data>, allocate memory and return pointer.
 */
M8_RESULTS* 
M8_RESULTS_Create();

/*! FUNCTION:  M8_RESULTS_Destroy()
 *  SYNOPSIS:  Destroy <results>, free memory and return NULL pointer.
 */
M8_RESULTS* 
M8_RESULTS_Destroy( M8_RESULTS* results );

/*! FUNCTION:  M8_RESULTS_Pushback()
 *  SYNOPSIS:  Add <res> to <results> list, resize array if full.
 */
void 
M8_RESULTS_Pushback(    M8_RESULTS*    results, 
                        M8_RESULT*     res );

/*! FUNCTION:  M8_RESULTS_Pushback()
 *  SYNOPSIS:  Resize <results> list to be able to store <size> entries.
 */
void 
M8_RESULTS_Resize(   M8_RESULTS* 	res,
                     size_t   	   size );

/*! FUNCTION:  M8_RESULTS_GetX()
 *  SYNOPSIS:  Get <i>th entry in <results> list.
 */
M8_RESULT* 
M8_RESULTS_GetX(  M8_RESULTS* 	results,
                  int   	      i );

/*! FUNCTION:  M8_RESULTS_Swap_Target_and_Query()
 *  SYNOPSIS:  The target and query are cross-labeled between MMSEQS and MMORE.
 *             Ideally, this should be remedied and MMORE should be swapped (todo list).
 *             For now, this will swap the fields so that m8 results are corrected for MMORE.
 */
M8_RESULT* 
M8_RESULTS_Swap_Target_and_Query(   M8_RESULTS* 	results );

/*! FUNCTION:  M8_RESULTS_Dump()
 *  SYNOPSIS:  Output all entries in <results> in .m8 format to file <fp>.
 */
void 
M8_RESULTS_Dump(  M8_RESULTS*    results,
                  FILE*          fp );

/*! FUNCTION:  M8_RESULTS_Dump()
 *  SYNOPSIS:  Output single <result> to file <fp>.
 */
void 
M8_RESULT_Dump(   M8_RESULT*  result,
                  FILE*       fp );

#endif /* _M8_RESULTS_H */