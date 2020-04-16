/*******************************************************************************
 *  @file cloud_search.h
 *  @brief Testing for navigating through the matrices.
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/


#ifndef _TESTING_H
#define _TESTING_H


/* test to cycle through all diags */
void fwd_test_cycle( const int   Q, 
                     const int   T,
                     MATRIX_3D*  st_MX,
                     MATRIX_2D*  sp_MX,
                     ALIGNMENT*  tr );

/* cycle through all indices in linear matrix, diag-by-diag */
void fwd_test_cycle3(const int   Q, 
                     const int   T, 
                     MATRIX_3D*  st_MX, 
                     MATRIX_3D*  st_MX3,
                     MATRIX_2D*  sp_MX, 
                     ALIGNMENT*  tr );

/* cycle through all indices in quadratic matrix in reverse, diag-by-diag */
void bck_test_cycle( const int   Q, 
                     const int   T,
                     MATRIX_3D*  st_MX,
                     MATRIX_2D*  sp_MX,
                     ALIGNMENT*  tr );

/* cycle through all indices in linear matrix in reverse, diag-by-diag */
void bck_test_cycle3(const int   Q, 
                     const int   T, 
                     MATRIX_3D*  st_MX, 
                     MATRIX_3D*  st_MX3,
                     MATRIX_2D*  sp_MX, 
                     ALIGNMENT*  tr );

/* test to show the cloud area, fill with value, return number of used cells  */
int cloud_Fill(const int   Q, 
               const int   T,
               MATRIX_3D*  st_MX,
               MATRIX_2D*  sp_MX,
               EDGEBOUNDS* edg,
               float       val, 
               int         mode );

/* test to show the cloud area, fill with value, return number of used cells  */
int cloud_Solid_Fill(const int      Q, 
                     const int      T,
                     MATRIX_3D*     st_MX,
                     MATRIX_2D*     sp_MX,
                     EDGEBOUNDS*    edg,
                     float          val, 
                     int            mode );

/* test to show the cloud area, fill with value, return number of used cells  */
int cloud_Cell_Count(const int   Q, 
                     const int   T,
                     MATRIX_3D*  st_MX,
                     MATRIX_2D*  sp_MX );

/* DP MATRIX FUNCTIONS */
void dp_matrix_Print(const int         Q, 
                     const int         T,
                     MATRIX_3D*        st_MX,
                     MATRIX_2D*        sp_MX );

void dp_matrix_Print3(  const int         T, 
                        const int         Q,
                        const MATRIX_3D*  st_MX3 );

void test_matrix_Print( const int         Q, 
                        const int         T,
                        const MATRIX_3D*  test_MX );

/* Clear all matrix values to -INF. (for testing) */
void dp_matrix_Clear(const int   Q, 
                     const int   T, 
                     MATRIX_3D*  st_MX, 
                     MATRIX_2D*  sp_MX );

/* Clear all matrix values to -INF. (for testing) */
void dp_matrix_Clear3(  const int   Q, 
                        const int   T,
                        MATRIX_3D*  st_MX3,
                        MATRIX_2D*  sp_MX );

/* Set all matrix values to val */
void dp_matrix_Clear_X( const int   Q, 
                        const int   T, 
                        MATRIX_3D*  st_MX, 
                        MATRIX_2D*  sp_MX,
                        float       val );

/* Clear all matrix values to -INF. (for testing) */
void dp_matrix_Clear_X3(const int   Q, 
                        const int   T,
                        MATRIX_3D*  st_MX3,
                        MATRIX_2D*  sp_MX,
                        int         val);

void dp_matrix_Save( const int         Q, 
                     const int         T, 
                     const MATRIX_3D*  st_MX, 
                     const MATRIX_2D*  sp_MX,
                     const char*       _filename_ );

void dp_matrix_trace_Save( const int         Q, 
                           const int         T, 
                           const MATRIX_3D*  st_MX, 
                           const MATRIX_2D*  sp_MX,
                           const ALIGNMENT*  tr,
                           const char*       _filename_ );

/* Copy source matrix into destination matrix */
void dp_matrix_Copy (const int   Q, 
                     const int   T,
                     MATRIX_3D*  st_MX_src,
                     MATRIX_2D*  sp_MX_src,
                     MATRIX_3D*  st_MX_dst,
                     MATRIX_2D*  sp_MX_dst );

int dp_matrix_Compare ( const int   Q, 
                        const int   T,
                        MATRIX_3D*  st_MX_1,
                        MATRIX_2D*  sp_MX_1,
                        MATRIX_3D*  st_MX_2,
                        MATRIX_2D*  sp_MX_2 );

void dp_matrix_Dump( const int         Q, 
                     const int         T, 
                     const MATRIX_3D*  st_MX, 
                     const MATRIX_2D*  sp_MX,
                     FILE*             fp );

void dp_matrix_trace_Save( const int         Q, 
                           const int         T, 
                           const MATRIX_3D*  st_MX, 
                           const MATRIX_2D*  sp_MX,
                           const ALIGNMENT*  tr,
                           const char*       _filename_ );

#endif /* _TESTING_H */