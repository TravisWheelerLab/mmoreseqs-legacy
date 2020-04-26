/*******************************************************************************
 *  FILE:      matrix_2d.h
 *  PURPOSE:   MATRIX_2D Float object.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _MATRIX_2D_H
#define _MATRIX_2D_H

/* constructor */
MATRIX_2D* MATRIX_2D_Create(int  R, 
                            int  C );

/* destructor */
void MATRIX_2D_Destroy(MATRIX_2D*  mx);

/* deep copy: returns dest matrix, will allocate if null */
MATRIX_2D* MATRIX_2D_Copy( MATRIX_2D*     dest,
                           MATRIX_2D*     src );

/* fill */
MATRIX_2D* MATRIX_2D_Fill(MATRIX_2D*  mx,
                          float       val);

/* getter for index */
float* MATRIX_2D_Get(MATRIX_2D*  mx, 
                     int         i, 
                     int         j);

/* convert 2D-coords to 1D-coords */
int MATRIX_2D_to_1D(MATRIX_2D*  mx, 
                    int         i, 
                    int         j);

/* reuse matrix by resizing only if new matrix requires more memory */
float MATRIX_2D_Reuse(MATRIX_2D*  mx, 
                      int         R, 
                      int         C);

/* resize matrix  */
float MATRIX_2D_Resize(MATRIX_2D*  mx, 
                       int         R, 
                       int         C);

/* Outputs MATRIX_2D out to FILE POINTER */
void MATRIX_2D_Dump(MATRIX_2D*  mx,
                    FILE*       fp);

/* Save MATRIX_2D to FILE by FILENAME */
void MATRIX_2D_Save( MATRIX_2D*  mx,
                     char*       _filename_ );

/* Compare two MATRIX_2D */
int MATRIX_2D_Compare(  MATRIX_2D*  mx_a,
                        MATRIX_2D*  mx_b );

/*
 *  FUNCTION:  MATRIX_2D_Diff()
 *  SYNOPSIS:  Takes difference of <mx_a> - <mx_b>.  Result stored in <mx_diff>.
 */
int MATRIX_2D_Diff(  MATRIX_2D*  mx_a,
                     MATRIX_2D*  mx_b,
                     MATRIX_2D*  mx_diff );

/* Unit Test */
void MATRIX_2D_UnitTest();

#endif /* _MATRIX_2D_H */