/*******************************************************************************
 *  FILE:      testing.h
 *  PURPOSE:   Testing for navigating through the matrices.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _TESTING_H
#define _TESTING_H

/*! FUNCTION:  TEST_set_color()
 *  SYNOPSIS:  Test output of all colors.
 */
void 
TEST_cycle_colors();

/*! FUNCTION:  TEST_set_color()
 *  SYNOPSIS:  Set console text color.
 */
void 
TEST_set_color(   char*    color,
                  bool     bold );

/*! FUNCTION:  TEST_set_color_num()
 *  SYNOPSIS:  Set console text color by index number <color> and boolean <bold>.
 */
void 
TEST_set_color_num(  int    color,
                     bool    bold );

/*! FUNCTION:  TEST_fwd_cycle()
 *  SYNOPSIS:  Cycle through all indices in quadratic matrix in forward direction, antidiag-by-antidiag 
 */
void 
TEST_fwd_cycle(   const int   Q, 
                  const int   T,
                  MATRIX_3D*  st_MX,
                  MATRIX_2D*  sp_MX,
                  ALIGNMENT*  tr );

/*! FUNCTION:  TEST_fwd_cycle3()
 *  SYNOPSIS:  Cycle through all indices in linear matrix, antidiag-by-antidiag 
 */
void 
TEST_fwd_cycle3(  const int   Q, 
                  const int   T, 
                  MATRIX_3D*  st_MX, 
                  MATRIX_3D*  st_MX3,
                  MATRIX_2D*  sp_MX, 
                  ALIGNMENT*  tr );

/*! FUNCTION:  TEST_bck_cycle()
 *  SYNOPSIS:  Cycle through all indices in quadratic matrix in backward direction, antidiag-by-antidiag 
 */
void 
TEST_bck_cycle(   const int   Q, 
                  const int   T,
                  MATRIX_3D*  st_MX,
                  MATRIX_2D*  sp_MX,
                  ALIGNMENT*  tr );

/*! FUNCTION:  TEST_bck_cycle3()
 *  SYNOPSIS:  Cycle through all indices in linear matrix in backward direction, antidiag-by-antidiag 
 */
void 
TEST_bck_cycle3(  const int   Q, 
                  const int   T, 
                  MATRIX_3D*  st_MX, 
                  MATRIX_3D*  st_MX3,
                  MATRIX_2D*  sp_MX, 
                  ALIGNMENT*  tr );

/*! FUNCTION:  MATRIX_2D_Cloud_Fill()
 *  SYNOPSIS:  Fill MATRIX_2D with value according to EDGEBOUNDS, returns number of cells covered by EDGEBOUNDS 
 */
int 
MATRIX_2D_Cloud_Fill(   MATRIX_2D*     cloud_MX,
                        EDGEBOUNDS*    edg,
                        float          val );

/*! FUNCTION:  MATRIX_2D_Cloud_Compare()
 *  SYNOPSIS:  Compare <cloud_MX_a> and <cloud_MX_b> ( equality is if cells are both positive, negative, or zero ).
 *             When debugging, stores heatmap of differences into <debugger->cloud_MX>.
 */
int 
MATRIX_2D_Cloud_Compare(   MATRIX_2D*     cloud_MX_a,
                           MATRIX_2D*     cloud_MX_b );

/*! FUNCTION:  MATRIX_2D_Cloud_Count()
 *  SYNOPSIS:  Count number of cells in MATRIX_2D with positive values. 
 */
int 
MATRIX_2D_Cloud_Count(  MATRIX_2D*  cloud_MX );

/*! FUNCTION: EDGEBOUNDS_Compare_by_Cloud()
 *  SYNOPSIS: Compare two EDGEBOUNDS by filling cloud matrices.
 */
int 
EDGEBOUNDS_Compare_by_Cloud(  EDGEBOUNDS*    edg_a,
                              MATRIX_2D*     mx_a,
                              EDGEBOUNDS*    edg_b,
                              MATRIX_2D*     mx_b );

/*! FUNCTION:  EDGEBOUNDS_Compare_by_Cloud_Single()
 *  SYNOPSIS:  Compares two EDGEBOUNDS by filling single cloud matrix.
 *             If equal, returns 0.  Else number of inequal cells.
 */
int 
EDGEBOUNDS_Compare_by_Cloud_Single(    MATRIX_2D*     mx,
                                       EDGEBOUNDS*    edg_a,
                                       EDGEBOUNDS*    edg_b );

#endif /* _TESTING_H */