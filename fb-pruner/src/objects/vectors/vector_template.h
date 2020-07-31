/*******************************************************************************
 *  FILE:      vector_template.c
 *  PURPOSE:   VECTOR_TMP Object Functions
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _VECTOR_TMP_H
#define _VECTOR_TMP_H

/*
 *  FUNCTION:  VECTOR_TMP_Create()
 *  SYNOPSIS:  Allocates new VECTOR_TMP object and returns pointer.
 */
VECTOR_TMP* VECTOR_TMP_Create();

/*
 *  FUNCTION:  VECTOR_TMP_Create()
 *  SYNOPSIS:  Allocates new VECTOR_TMP object at specific size and returns pointer.
 */
VECTOR_TMP* VECTOR_TMP_Create_by_Size( int    size );

/*
 *  FUNCTION:  VECTOR_TMP_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_TMP.
 */
void* VECTOR_TMP_Destroy( VECTOR_TMP*   vec );

/*
 *  FUNCTION:  VECTOR_TMP_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_TMP object by resetting size counter (no realloc) .
 */
void VECTOR_TMP_Reuse( VECTOR_TMP*   vec );

/*
 *  FUNCTION:  VECTOR_TMP_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object. 
 *             Creates new VECTOR_TMP for <dest> if <dest> is NULL.
 */
VECTOR_TMP* VECTOR_TMP_Copy(  VECTOR_TMP*   src, 
                              VECTOR_TMP*   dest );

/*
 *  FUNCTION:  VECTOR_TMP_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>. 
 */
void VECTOR_TMP_Resize(    VECTOR_TMP*    vec, 
                           int            size );

/*
 *  FUNCTION:  VECTOR_TMP_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>. 
 */
void VECTOR_TMP_GrowTo(   VECTOR_TMP*     vec, 
                           int            size );

/*
 *  FUNCTION:  VECTOR_TMP_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
void VECTOR_TMP_Pushback(  VECTOR_TMP*   vec, 
                           TMP           val );

/*
 *  FUNCTION:  VECTOR_TMP_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data. 
 */
TMP VECTOR_TMP_Pop( VECTOR_TMP*   vec );

/*
 *  FUNCTION:  VECTOR_TMP_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
void VECTOR_TMP_Set(   VECTOR_TMP*     vec, 
                        int            idx, 
                        TMP            val );

/*
 *  FUNCTION:  VECTOR_TMP_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
TMP VECTOR_TMP_Get(  VECTOR_TMP*   vec, 
                     int           idx );

/*
 *  FUNCTION:  VECTOR_TMP_Get_Ref()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return pointer to data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
TMP* VECTOR_TMP_Get_Ref(   VECTOR_TMP*   vec, 
                           int           idx );

/*
 *  FUNCTION:  VECTOR_TMP_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality, 
 *             pos if <vec_A> > <vec_B>,  
 *             neg if <vec_A> < <vec_B>.
 */
int VECTOR_TMP_Compare(    VECTOR_TMP*   vec_A, 
                           VECTOR_TMP*   vec_B );

/*
 *  FUNCTION:  VECTOR_TMP_Sort()
 *  SYNOPSIS:  Sort <vec> data array in ascending order.
 */
void VECTOR_TMP_Sort( VECTOR_TMP*    vec );

/*
 *  FUNCTION:  VECTOR_TMP_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
void VECTOR_TMP_Dump(  VECTOR_TMP*   vec,
                           FILE*             fp );

#endif 