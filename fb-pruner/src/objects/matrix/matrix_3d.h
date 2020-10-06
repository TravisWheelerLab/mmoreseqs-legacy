/*******************************************************************************
 *  FILE:      matrix_3d.c
 *  PURPOSE:   MATRIX_3D Float object.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _MATRIX_3D_H
#define _MATRIX_3D_H

/** FUNCTION:  MATRIX_3D_Create()
 *  SYNOPSIS:  
 */
MATRIX_3D* 
MATRIX_3D_Create( const int  R,
                  const int  C,
                  const int  N );


/** FUNCTION:  MATRIX_3D_Destroy()
 *  SYNOPSIS:  
 */
MATRIX_3D* 
MATRIX_3D_Destroy( MATRIX_3D*  mx );


/** FUNCTION:  MATRIX_3D_Create_Clean()
 *  SYNOPSIS:  
 */
MATRIX_3D* 
MATRIX_3D_Create_Clean( const int  R,
                        const int  C,
                        const int  N );


/** FUNCTION:  MATRIX_3D_Copy()
 *  SYNOPSIS:  
 */
MATRIX_3D* 
MATRIX_3D_Copy(   MATRIX_3D*           dest,
                  const MATRIX_3D*     src );


/** FUNCTION:  MATRIX_3D_Fill()
 *  SYNOPSIS:  
 */
void 
MATRIX_3D_Fill(   MATRIX_3D*     mx,
                  const float    val );


/** FUNCTION:  MATRIX_3D_Clean()
 *  SYNOPSIS:  
 */
void 
MATRIX_3D_Clean(   MATRIX_3D*     mx );


/** FUNCTION:  MATRIX_3D_Check_Clean()
 *  SYNOPSIS:  
 */
int 
MATRIX_3D_Check_Clean( MATRIX_3D*     mx );


/** FUNCTION:  MATRIX_3D_Check_Value()
 *  SYNOPSIS:  
 */
int 
MATRIX_3D_Check_Value(  MATRIX_3D*     mx,
                        const float    val );


/** FUNCTION:  MATRIX_3D_Get()
 *  SYNOPSIS:  
 */
float* 
MATRIX_3D_Get( MATRIX_3D*  mx,
               const int   i,
               const int   j,
               const int   k );


/** FUNCTION:  MATRIX_3D_Get_X()
 *  SYNOPSIS:  Get pointer to cell from <mx> at position <i,j,k>.
 */
float* 
MATRIX_3D_Get_X(  MATRIX_3D*  mx,
                  const int   i,
                  const int   j,
                  const int   k );


/** FUNCTION:  MATRIX_3D_Get_1D()
 *  SYNOPSIS:  Get pointer to cell from <mx> at flat-index position <n>.
 */
float* 
MATRIX_3D_Get_1D( MATRIX_3D*  mx,
                  const int   n );


/** FUNCTION:  MATRIX_3D_to_1D()
 *  SYNOPSIS:  
 */
int 
MATRIX_3D_to_1D(  const MATRIX_3D*  mx,
                  const int         i,
                  const int         j,
                  const int         k );


/** FUNCTION:  MATRIX_3D_Reuse()
 *  SYNOPSIS:  Reuse MATRIX_3D by resizing only if new matrix requires more memory.
 */
float 
MATRIX_3D_Reuse(  MATRIX_3D*  mx,
                  const int   R,
                  const int   C,
                  const int   N );


/** FUNCTION:  MATRIX_3D_Reuse_Clean()
 *  SYNOPSIS:  Reuse MATRIX_3D by resizing only if new matrix requires more memory. 
 *             New memory is set to -INF.
 */
float 
MATRIX_3D_Reuse_Clean(  MATRIX_3D*  mx,
                        const int   R,
                        const int   C,
                        const int   N );


/** FUNCTION:  MATRIX_3D_Resize()
 *  SYNOPSIS:  Resize MATRIX_3D to new dimensions.
 */
float 
MATRIX_3D_Resize( MATRIX_3D*  mx,
                  const int   R,
                  const int   C,
                  const int   N );


/** FUNCTION:  MATRIX_3D_Dump()
 *  SYNOPSIS:  Outputs MATRIX_3D <mx> out to file pointer <fp>.
 */
void 
MATRIX_3D_Dump(   MATRIX_3D*  mx,
                  FILE*       fp );


/** FUNCTION:  MATRIX_3D_Save()
 *  SYNOPSIS:  Save MATRIX_3D <mx> to file with filename <filename>.
 */
void 
MATRIX_3D_Save(   MATRIX_3D*  mx,
                  char*       _filename_);


/** FUNCTION:  MATRIX_3D_Compare()
 *  SYNOPSIS:  Compare MATRIX_3D's <mx_A> and <mx_B>.
 */
int 
MATRIX_3D_Compare(   MATRIX_3D*  mx_a,
                     MATRIX_3D*  mx_b );

/** FUNCTION:  MATRIX_3D_Add()
 *  SYNOPSIS:  Takes sum of <mx_a> + <mx_b>.  Result stored in <mx_res>.
 */
void 
MATRIX_3D_Add( MATRIX_3D*  mx_a,
               MATRIX_3D*  mx_b,
               MATRIX_3D*  mx_res );

/** FUNCTION:  MATRIX_3D_Diff()
 *  SYNOPSIS:  Takes difference of <mx_a> - <mx_b>.  Result stored in <mx_diff>.
 */
void 
MATRIX_3D_Diff(   MATRIX_3D*  mx_a,
                  MATRIX_3D*  mx_b,
                  MATRIX_3D*  mx_diff );

/** FUNCTION:  MATRIX_3D_Utest()
 *  SYNOPSIS:  Unit Test for MATIX_3D.
 */
void MATRIX_3D_Utest();

#endif /* _MATRIX_3D_H */