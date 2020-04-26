/*******************************************************************************
 *  FILE:      matrix_3d.c
 *  PURPOSE:   MATRIX_3D Float object.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _MATRIX_3D_H
#define _MATRIX_3D_H

/* constructor */
MATRIX_3D* MATRIX_3D_Create(const int  R,
                            const int  C,
                            const int  N);
/* destructor */
void MATRIX_3D_Destroy(MATRIX_3D*  mx);

/* deep copy */
MATRIX_3D* MATRIX_3D_Copy( MATRIX_3D*     dest,
                           MATRIX_3D*     src );

/* getter for index */
float* MATRIX_3D_Get(MATRIX_3D*  mx,
                     const int   i,
                     const int   j,
                     const int   k );

/* convert 3D-coords to 1D-coords */
int MATRIX_3D_to_1D(const MATRIX_3D*  mx,
                    const int         i,
                    const int         j,
                    const int         k );

/* reuse MATRIX_3D by resizing only if new matrix requires more memory */
float MATRIX_3D_Reuse(MATRIX_3D*  mx,
                      const int   R,
                      const int   C,
                      const int   N );

/* fill MATRIX_3D with values */
void MATRIX_3D_Fill( MATRIX_3D*  mx,
                     const float val );

/* resize MATRIX_3D to new dimensions */
float MATRIX_3D_Resize(MATRIX_3D*  mx,
                       const int   R,
                       const int   C,
                       const int   N);

/* Outputs MATRIX_3D out to FILE POINTER */
void MATRIX_3D_Dump(MATRIX_3D*  mx,
                    FILE*       fp);

/* Save MATRIX_3D to FILE by FILENAME */
void MATRIX_3D_Save(MATRIX_3D*  mx,
                    char*       _filename_);

/* Compare two MATRIX_2D */
int MATRIX_3D_Compare(  MATRIX_3D*  mx_a,
                        MATRIX_3D*  mx_b );

/*
 *  FUNCTION:  MATRIX_2D_Diff()
 *  SYNOPSIS:  Takes difference of <mx_a> - <mx_b>.  Result stored in <mx_diff>.
 */
void MATRIX_3D_Diff( MATRIX_3D*  mx_a,
                     MATRIX_3D*  mx_b,
                     MATRIX_3D*  mx_diff );

/* unit test */
void MATRIX_3D_Utest();

#endif /* _MATRIX_3D_H */