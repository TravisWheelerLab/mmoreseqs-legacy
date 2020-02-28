/******************************************************************************
 *  @file matrix.c
 *  @brief float matrix object
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

#ifndef _MATRIX_2D_H
#define _MATRIX_2D_H

// typedef struct {
//    int    R; 
//    int    C;
//    int    Nalloc;
//    float* data;
// } MATRIX_2D;

/* constructor */
MATRIX_2D* MATRIX_2D_Create(int  R, 
                            int  C);

/* destructor */
void MATRIX_2D_Destroy(MATRIX_2D*  mx);

/* fill */
MATRIX_2D* MATRIX_2D_Fill(MATRIX_2D*  mx,
                          int         R, 
                          int         C, 
                          int         val);

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

/* Unit Test */
void MATRIX_2D_UnitTest();

#endif /* _MATRIX_2D_H */