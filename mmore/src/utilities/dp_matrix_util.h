/*******************************************************************************
 *  FILE:      dp_matrix.h
 *  SYNOPSIS:  Dynamic Programming Matrix functions.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _DP_MATRIX_UTIL_H
#define _DP_MATRIX_UTIL_H

/* === imports === */
#include "../objects/structs.h"

/* === macros === */
/* NONE */

/* === public functions === */

/*! FUNCTION:  DP_MATRIX_GetBounds()
 *  SYNOPSIS:  Get the edgebounds of matrix at given antidiagonal (closed form).
 *             Stores new bounds in BOUND <bnd>.
 */
void 
DP_MATRIX_GetBounds(   const int   Q,
                        const int   T,
                        int         d_0,
                        int         d_cnt,
                        BOUND*      bnd );

/*! FUNCTION:  DP_MATRIX_Copy()
 *  SYNOPSIS:  Copy dynamic programming matrix into destination.
 */
void 
DP_MATRIX_Copy(   const int      Q, 
                  const int      T,
                  MATRIX_3D*     st_MX_src,
                  MATRIX_2D*     sp_MX_src,
                  MATRIX_3D*     st_MX_dst,
                  MATRIX_2D*     sp_MX_dst );

/*! FUNCTION:  DP_MATRIX_Fill()
 *  SYNOPSIS:  Fill entire dynamic programming matrix with value
 */
void 
DP_MATRIX_Fill(   const int      Q, 
                  const int      T,
                  MATRIX_3D*     st_MX,
                  MATRIX_2D*     sp_MX,
                  const float    val );

/*! FUNCTION:  DP_MATRIX_Clean()
 *  SYNOPSIS:  Fill entire dynamic programming matrix with clean.
 */
void 
DP_MATRIX_Clean(  const int   Q, 
                  const int   T,
                  MATRIX_3D*  st_MX,
                  MATRIX_2D*  sp_MX );

/*! FUNCTION:  DP_MATRIX_Clean_Verify()
 *  SYNOPSIS:  Check whether there are clean.  If clean, returns true.
 */
bool 
DP_MATRIX_Clean_Verify(    const int   Q, 
                           const int   T,
                           MATRIX_3D*  st_MX,
                           MATRIX_2D*  sp_MX );

/*! FUNCTION:  DP_MATRIX_Save()
 *  SYNOPSIS: Compare two dynamic programming matrices.
 */
int 
DP_MATRIX_Compare(   MATRIX_3D*  st_MX_1,
                     MATRIX_2D*  sp_MX_1,
                     MATRIX_3D*  st_MX_2,
                     MATRIX_2D*  sp_MX_2 );


/*! FUNCTION:  DP_MATRIX_Diff()
 *  SYNOPSIS: Compare two dynamic programming matrices.
 */
int 
DP_MATRIX_Diff(   MATRIX_3D*  st_MX_1,
                  MATRIX_2D*  sp_MX_1,
                  MATRIX_3D*  st_MX_2,
                  MATRIX_2D*  sp_MX_2,
                  MATRIX_3D*  st_MX_res,
                  MATRIX_2D*  sp_MX_res );

/*! FUNCTION:  DP_MATRIX_Save()
 *  SYNOPSIS:  Save dynamic programming matrix to file (by filename).
 */
void 
DP_MATRIX_Save(   const int         Q, 
                  const int         T, 
                  MATRIX_3D*        st_MX, 
                  MATRIX_2D*        sp_MX,
                  const char*       filename );

/*! FUNCTION:  DP_MATRIX_Dump()
 *  SYNOPSIS:  Output dynamic programming matrix to file.
 */
void 
DP_MATRIX_Dump(   const int         Q, 
                  const int         T, 
                  MATRIX_3D*        st_MX, 
                  MATRIX_2D*        sp_MX,
                  FILE*             fp );

/*! FUNCTION:  DP_MATRIX_Log_Dump()
 *  SYNOPSIS:  Output dynamic programming matrix to file.
 */
void 
DP_MATRIX_Log_Dump(  const int         Q,
                     const int         T,
                     MATRIX_3D*        st_MX,
                     MATRIX_2D*        sp_MX,
                     FILE*             fp );

/*! FUNCTION:  DP_MATRIX_Dump()
 *  SYNOPSIS:  Output dynamic programming matrix to file.
 */
void 
DP_MATRIX_Norm_Dump(    const int         Q, 
                        const int         T, 
                        MATRIX_3D*        st_MX, 
                        MATRIX_2D*        sp_MX,
                        FILE*             fp );

/*! FUNCTION:  DP_MATRIX_Dump()
 *  SYNOPSIS:  Output dynamic programming matrix to file.
 */
void DP_MATRIX_Spec_Dump(  const int         Q, 
                           const int         T, 
                           MATRIX_3D*        st_MX, 
                           MATRIX_2D*        sp_MX,
                           FILE*             fp );

/*! FUNCTION:  DP_MATRIX_Dump()
 *  SYNOPSIS:  Output dynamic programming matrix to file.
 */
void 
DP_MATRIX_SpecExp_Dump(    const int         Q,
                           const int         T,
                           MATRIX_3D*        st_MX,
                           MATRIX_2D*        sp_MX,
                           FILE*             fp );

/*! FUNCTION:  DP_MATRIX_Trace_Save()
 *  SYNOPSIS:  Save dynamic programming matrix with viterbi trace to file.
 */
void 
DP_MATRIX_Trace_Save(   const int         Q, 
                        const int         T, 
                        MATRIX_3D*        st_MX, 
                        MATRIX_2D*        sp_MX,
                        ALIGNMENT*        tr,
                        const char*       filename );

/*! FUNCTION:  DP_MATRIX_Trace_Save()
 *  SYNOPSIS:  Output dynamic programming matrix with viterbi trace to file.
 */
void 
DP_MATRIX_Trace_Dump(   const int         Q, 
                        const int         T, 
                        MATRIX_3D*        st_MX, 
                        MATRIX_2D*        sp_MX,
                        ALIGNMENT*        tr,
                        FILE*             fp );

/*! FUNCTION:  DP_MATRIX_VIZ_Compare()
 *  SYNOPSIS:  Projects two EDGEBOUNDS onto 2D_MATRIX.
 */
void 
DP_MATRIX_VIZ_Compare(  MATRIX_2D*        cloud_MX,
                        EDGEBOUNDS*       edg_1,
                        EDGEBOUNDS*       edg_2 );

/*! FUNCTION:  DP_MATRIX_VIZ_Save()
 *  SYNOPSIS:  Saves simple visualization to filename.
 */
void 
DP_MATRIX_VIZ_Save(  MATRIX_2D*     cloud_MX,
                     char*          filename );

/*! FUNCTION:  DP_MATRIX_VIZ_Dump()
 *  SYNOPSIS:  Outputs simple visualization of dp matrix.
 */
void 
DP_MATRIX_VIZ_Dump(  MATRIX_2D*        test_MX,
                     FILE*             fp );

/*! FUNCTION:  DP_MATRIX_VIZ_Dump()
 *  SYNOPSIS:  Outputs simple visualization of dp matrix.
 *             Symbols are color-coded.
 */
void 
DP_MATRIX_VIZ_Color_Dump(  MATRIX_2D*        cloud_MX,
                           FILE*             fp );

/*! FUNCTION:  DP_MATRIX_VIZ_Trace()
 *  SYNOPSIS:  Adds trace to visulization of dp matrix.
 */
void 
DP_MATRIX_VIZ_Trace(    MATRIX_2D*        cloud_MX,
                        const ALIGNMENT*  aln );

/*! FUNCTION:  DP_MATRIX_MAT_Dump()
 *  SYNOPSIS:  Outputs match state of matrix.
 */
void 
DP_MATRIX_MAT_Dump(  int         Q,
                     int         T,
                     MATRIX_3D*  st_MX,
                     FILE*       fp );

/*! FUNCTION: EDGEBOUNDS_BruteCount()
 *  SYNOPSIS: For testing purposes to get the true count of edgebounds. 
 *            Fills 2D_MATRIX will zeros, then fills with edgebounds with ones.
 *            Sums over matrix.
 */
int
EDGEBOUNDS_BruteCount( EDGEBOUNDS* edg );

#endif /* _DP_MATRIX_UTIL_H */