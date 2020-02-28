/*******************************************************************************
 *  @file matrix.c
 *  @brief float matrix object
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

#ifndef _MATRIX_3D_H
#define _MATRIX_3D_H

// typedef struct {
//    int     R; 
//    int     C;
//    int     N;
//    int     Nalloc;
//    float*  data;
// } MATRIX_3D;

/* constructor */
MATRIX_3D* MATRIX_3D_Create(int  R, 
                            int  C,
                            int  N);
/* destructor */
void MATRIX_3D_Destroy(MATRIX_3D*  mx);

/* fill */
MATRIX_3D* MATRIX_3D_Fill(MATRIX_3D*  mx,
                          int         R, 
                          int         C,
                          int         N, 
                          int         val);
/* getter for index */
float* MATRIX_3D_Get(MATRIX_3D*  mx, 
                     int         i, 
                     int         j,
                     int         k);
/* convert 3D-coords to 1D-coords */
int MATRIX_3D_to_1D(MATRIX_3D*  mx, 
                    int         i, 
                    int         j,
                    int         k);
/* reuse matrix by resizing only if new matrix requires more memory */
float MATRIX_3D_Reuse(MATRIX_3D*  mx, 
                      int         R, 
                      int         C,
                      int         N);
/* resize MATRIX_3D to new dimensions */
float MATRIX_3D_Resize(MATRIX_3D*  mx, 
                       int         R, 
                       int         C,
                       int         N);
/* Outputs MATRIX_3D out to FILE POINTER */
void MATRIX_3D_Dump(MATRIX_3D*  mx,
                    FILE*       fp);

/* unit test */
void MATRIX_3D_UnitTest();

#endif /* _MATRIX_3D_H */