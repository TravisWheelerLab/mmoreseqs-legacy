/*******************************************************************************
 *  FILE:      vector_char.c
 *  PURPOSE:   VECTOR_CHAR Object Functions
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _VECTOR_CHAR_H
#define _VECTOR_CHAR_H

/*
 *  FUNCTION:  VECTOR_CHAR_Create()
 *  SYNOPSIS:  Allocates new VECTOR_CHAR object and returns pointer.
 */
VECTOR_CHAR* VECTOR_CHAR_Create();

/*
 *  FUNCTION:  VECTOR_CHAR_Create()
 *  SYNOPSIS:  Allocates new VECTOR_CHAR object at specific size and returns pointer.
 */
VECTOR_CHAR* VECTOR_CHAR_Create_by_Size( int    size );

/*
 *  FUNCTION:  VECTOR_CHAR_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_CHAR.
 */
void* VECTOR_CHAR_Destroy( VECTOR_CHAR*   vec );

/*
 *  FUNCTION:  VECTOR_CHAR_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_CHAR object by resetting size counter (no realloc) .
 */
void VECTOR_CHAR_Reuse( VECTOR_CHAR*   vec );

/*
 *  FUNCTION:  VECTOR_CHAR_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object. 
 *             Creates new VECTOR_CHAR for <dest> if <dest> is NULL.
 */
VECTOR_CHAR* VECTOR_CHAR_Copy(  VECTOR_CHAR*   src, 
                              VECTOR_CHAR*   dest );

/*
 *  FUNCTION:  VECTOR_CHAR_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>. 
 */
void VECTOR_CHAR_Resize(    VECTOR_CHAR*    vec, 
                           int            size );

/*
 *  FUNCTION:  VECTOR_CHAR_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>. 
 */
void VECTOR_CHAR_GrowTo(   VECTOR_CHAR*     vec, 
                           int            size );

/*
 *  FUNCTION:  VECTOR_CHAR_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
void VECTOR_CHAR_Pushback(  VECTOR_CHAR*   vec, 
                           CHAR           val );

/*
 *  FUNCTION:  VECTOR_CHAR_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data. 
 */
CHAR VECTOR_CHAR_Pop( VECTOR_CHAR*   vec );

/*
 *  FUNCTION:  VECTOR_CHAR_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
void VECTOR_CHAR_Set(   VECTOR_CHAR*     vec, 
                        int            idx, 
                        CHAR            val );

/*
 *  FUNCTION:  VECTOR_CHAR_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
CHAR VECTOR_CHAR_Get(  VECTOR_CHAR*   vec, 
                     int           idx );

/*
 *  FUNCTION:  VECTOR_CHAR_Get_Ref()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return pointer to data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
CHAR* VECTOR_CHAR_Get_Ref(   VECTOR_CHAR*   vec, 
                           int           idx );

/*
 *  FUNCTION:  VECTOR_CHAR_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality, 
 *             pos if <vec_A> > <vec_B>,  
 *             neg if <vec_A> < <vec_B>.
 */
int VECTOR_CHAR_Compare(    VECTOR_CHAR*   vec_A, 
                           VECTOR_CHAR*   vec_B );

/*
 *  FUNCTION:  VECTOR_CHAR_Sort()
 *  SYNOPSIS:  Sort <vec> data array in ascending order.
 */
void VECTOR_CHAR_Sort( VECTOR_CHAR*    vec );

/*
 *  FUNCTION:  VECTOR_CHAR_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
void VECTOR_CHAR_Dump(  VECTOR_CHAR*   vec,
                           FILE*             fp );

#endif 