/*******************************************************************************
 *  @file cloud_search.h
 *  @brief Testing for navigating through the matrices.
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/


#ifndef _TESTING_H
#define _TESTING_H

/* === INCLUDES === */
// #include "objects/structs.h"
// #include "objects/edgebound.h"
// #include "objects/alignment.h"


/* === FUNCTIONS === */
void fwd_test_cycle( const int   Q, 
                     const int   T,
                     float*      st_MX,
                     float*      sp_MX,
                     ALIGNMENT*  tr );

void bck_test_cycle( const int   Q, 
                     const int   T,
                     float*      st_MX,
                     float*      sp_MX,
                     ALIGNMENT*  tr );

void fwd_test_cycle3(const int   Q, 
                     const int   T, 
                     float*      st_MX, 
                     float*      st_MX3,
                     float*      sp_MX, 
                     ALIGNMENT*  tr );

void bck_test_cycle3(const int   Q, 
                     const int   T, 
                     float*      st_MX, 
                     float*      st_MX3,
                     float*      sp_MX, 
                     ALIGNMENT*  tr );

int cloud_Fill(const int   Q, 
               const int   T,
               float*      st_MX,
               float*      sp_MX,
               EDGEBOUNDS* edg,
               float       val, 
               int         mode );

int cloud_Solid_Fill(const int      Q, 
                     const int      T,
                     float*         st_MX,
                     float*         sp_MX,
                     EDGEBOUNDS*    edg,
                     float          val, 
                     int            mode );

int cloud_Cell_Count(const int   Q, 
                     const int   T,
                     float*      st_MX,
                     float*      sp_MX );

/* DP MATRIX FUNCTIONS */
void dp_matrix_Print(const int      Q, 
                     const int      T, 
                     const float*   st_MX, 
                     const float*   sp_MX );
void dp_matrix_Print3(const int     Q, 
                      const int     T, 
                      const float*  st_MX3 );
void dp_matrix_Clear(const int   Q, 
                     const int   T, 
                     float*      st_MX, 
                     float*      sp_MX );
void dp_matrix_Clear3(  const int   Q, 
                        const int   T,
                        float*      st_MX3,
                        float*      sp_MX );
void dp_matrix_Clear_X( const int   Q, 
                        const int   T, 
                        float*      st_MX, 
                        float*      sp_MX,
                        float       val );
void dp_matrix_Clear_X3(const int   Q, 
                        const int   T,
                        float*      st_MX3,
                        float*      sp_MX,
                        int         val);
void dp_matrix_Save( const int      Q, 
                     const int      T, 
                     const float*   st_MX, 
                     const float*   sp_MX,
                     const char*    _filename_ );
void dp_matrix_trace_Save( const int         Q, 
                           const int         T, 
                           const float*      st_MX, 
                           const float*      sp_MX,
                           const ALIGNMENT*  tr,
                           const char*       _filename_ );
void test_matrix_Print( const int      Q, 
                        const int      T, 
                        const float*   test_MX );
int dp_matrix_Compare(  const int      Q, 
                        const int      T,
                        float*         st_MX_1,
                        float*         sp_MX_1,
                        float*         st_MX_2,
                        float*         sp_MX_2 );
void dp_matrix_Copy( const int   Q, 
                     const int   T,
                     float*      st_MX_src,
                     float*      sp_MX_src,
                     float*      st_MX_dst,
                     float*      sp_MX_dst );

#endif /* _TESTING_H */