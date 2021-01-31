/*******************************************************************************
 *  FILE:      dp_matrix.c
 *  PURPOSE:   DP_MATRIX object.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _DP_MATRIX_H
#define _DP_MATRIX_H

/** FUNCTION:  DP_MATRIX_Create()
 *  SYNOPSIS:  Create new DP_MATRIX object and returns pointer.
 *             Most data is left NULL to be supplied by WORK_init().
 */
DP_MATRIX* 
DP_MATRIX_Create( DPMX_MODE mode );

/** FUNCTION:  DP_MATRIX_Create()
 *  SYNOPSIS:  Create new DP_MATRIX object and returns pointer.
 *             Most data is left NULL to be supplied by WORK_init().
 */
DP_MATRIX* 
DP_MATRIX_GrowTo(    DP_MATRIX*    dom_def,
                     int            size );

/** FUNCTION:  DP_MATRIX_Destroy()
 *  SYNOPSIS:  Frees DP_MATRIX object and returns pointer.
 */
DP_MATRIX* 
DP_MATRIX_Destroy( DP_MATRIX* dom_def );

#endif /* _DP_MATRIX_H */