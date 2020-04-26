/*******************************************************************************
 *  FILE:      dp_matrix.h
 *  SYNOPSIS:  Dynamic Programming Matrix functions.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _DP_MATRIX_H
#define _DP_MATRIX_H

/*
 *  FUNCTION:  DP_MATRIX_Get_Bounds()
 *  SYNOPSIS:  Get the edgebounds of matrix at given antidiagonal (closed form).
 *             Stores new bounds in BOUND <bnd>.
 */
void DP_MATRIX_Get_Bounds( const int   Q,
                           const int   T,
                           int         d_0,
                           int         d_cnt,
                           BOUND*      bnd );

/*
 *  FUNCTION:  DP_MATRIX_Copy()
 *  SYNOPSIS:  Copy dynamic programming matrix into destination.
 */
void DP_MATRIX_Copy (const int   Q, 
                     const int   T,
                     MATRIX_3D*  st_MX_src,
                     MATRIX_2D*  sp_MX_src,
                     MATRIX_3D*  st_MX_dst,
                     MATRIX_2D*  sp_MX_dst );

/*
 *  FUNCTION:  DP_MATRIX_Fill()
 *  SYNOPSIS:  Fill entire dynamic programming matrix with value
 */
void DP_MATRIX_Fill (const int   Q, 
                     const int   T,
                     MATRIX_3D*  st_MX,
                     MATRIX_2D*  sp_MX,
                     const float val );

/*
 *  FUNCTION:  DP_MATRIX_Save()
 *  SYNOPSIS: Compare two dynamic programming matrices.
 */
int DP_MATRIX_Compare ( MATRIX_3D*  st_MX_1,
                        MATRIX_2D*  sp_MX_1,
                        MATRIX_3D*  st_MX_2,
                        MATRIX_2D*  sp_MX_2 );


/*
 *  FUNCTION:  DP_MATRIX_Diff()
 *  SYNOPSIS: Compare two dynamic programming matrices.
 */
int DP_MATRIX_Diff ( MATRIX_3D*  st_MX_1,
                     MATRIX_2D*  sp_MX_1,
                     MATRIX_3D*  st_MX_2,
                     MATRIX_2D*  sp_MX_2,
                     MATRIX_3D*  st_MX_res,
                     MATRIX_2D*  sp_MX_res );

/*
 *  FUNCTION:  DP_MATRIX_Save()
 *  SYNOPSIS:  Save dynamic programming matrix to file (by filename).
 */
void DP_MATRIX_Save( const int         Q, 
                     const int         T, 
                     MATRIX_3D*        st_MX, 
                     MATRIX_2D*        sp_MX,
                     const char*       _filename_ );

/*
 *  FUNCTION:  DP_MATRIX_Dump()
 *  SYNOPSIS:  Output dynamic programming matrix to file.
 */
void DP_MATRIX_Dump( const int         Q, 
                     const int         T, 
                     MATRIX_3D*        st_MX, 
                     MATRIX_2D*        sp_MX,
                     FILE*             fp );

/*
 *  FUNCTION:  DP_MATRIX_Trace_Save()
 *  SYNOPSIS:  Save dynamic programming matrix with viterbi trace to file.
 */
void DP_MATRIX_Trace_Save( const int         Q, 
                           const int         T, 
                           MATRIX_3D*        st_MX, 
                           MATRIX_2D*        sp_MX,
                           ALIGNMENT*        tr,
                           const char*       _filename_ );

/*
 *  FUNCTION:  DP_MATRIX_Trace_Save()
 *  SYNOPSIS:  Output dynamic programming matrix with viterbi trace to file.
 */
void DP_MATRIX_Trace_Dump( const int         Q, 
                           const int         T, 
                           MATRIX_3D*        st_MX, 
                           MATRIX_2D*        sp_MX,
                           ALIGNMENT*        tr,
                           FILE*             fp );

/*
 *  FUNCTION:  DP_MATRIX_VIZ_Dump()
 *  SYNOPSIS:  Outputs simple visualization of dp matrix.
 */
void DP_MATRIX_VIZ_Dump(   MATRIX_2D*        test_MX,
                           FILE*             fp );

/*
 *  FUNCTION:  DP_MATRIX_VIZ_Dump()
 *  SYNOPSIS:  Outputs simple visualization of dp matrix.
 *             Symbols are color-coded.
 */
void DP_MATRIX_VIZ_Color_Dump(   MATRIX_2D*        cloud_MX,
                                 FILE*             fp );

/*
 *  FUNCTION:  DP_MATRIX_VIZ_Trace()
 *  SYNOPSIS:  Adds trace to visulization of dp matrix.
 */
void DP_MATRIX_VIZ_Trace(  MATRIX_2D*        cloud_MX,
                           const ALIGNMENT*  aln );

/*
 *  FUNCTION:  DP_MATRIX_MAT_Dump()
 *  SYNOPSIS:  Outputs match state of matrix.
 */
void DP_MATRIX_MAT_Dump(   int         Q,
                           int         T,
                           MATRIX_3D*  st_MX,
                           FILE*       fp );

#endif /* _DP_MATRIX_H */