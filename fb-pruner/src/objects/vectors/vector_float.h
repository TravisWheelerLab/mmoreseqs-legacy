/*******************************************************************************
 *  FILE:      vector_int.c
 *  PURPOSE:   VECTOR_FLT Object Functions
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _VECTOR_FLT_H
#define _VECTOR_FLT_H

/*
 *  FUNCTION:  VECTOR_FLT_Create()
 *  SYNOPSIS:  Allocates new VECTOR_FLT object and returns pointer.
 */
VECTOR_FLT* VECTOR_FLT_Create();

/*
 *  FUNCTION:  VECTOR_FLT_Create()
 *  SYNOPSIS:  Allocates new VECTOR_FLT object at specific size and returns pointer.
 */
VECTOR_FLT* VECTOR_FLT_Create_by_Size( int    size );

/*
 *  FUNCTION:  VECTOR_FLT_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_FLT.
 */
void* VECTOR_FLT_Destroy( VECTOR_FLT*   vec );

/*
 *  FUNCTION:  VECTOR_FLT_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_FLT object by resetting size counter (no realloc) .
 */
void VECTOR_FLT_Reuse( VECTOR_FLT*   vec );

/*
 *  FUNCTION:  VECTOR_FLT_Fill()
 *  SYNOPSIS:  Fill VECTOR_FLT object with val.
 */
void VECTOR_FLT_Fill(   VECTOR_FLT*   vec, 
                        FLT           val );

/*
 *  FUNCTION:  VECTOR_FLT_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object. 
 *             Creates new VECTOR_FLT for <dest> if <dest> is NULL.
 */
VECTOR_FLT* VECTOR_FLT_Copy(  VECTOR_FLT*   src, 
                              VECTOR_FLT*   dest );

/*
 *  FUNCTION:  VECTOR_FLT_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>. 
 */
void VECTOR_FLT_Resize(    VECTOR_FLT*    vec, 
                           int            size );

/*
 *  FUNCTION:  VECTOR_FLT_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>. 
 */
void VECTOR_FLT_GrowTo(   VECTOR_FLT*     vec, 
                           int            size );

/*
 *  FUNCTION:  VECTOR_FLT_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
void VECTOR_FLT_Pushback(  VECTOR_FLT*   vec, 
                           FLT           val );

/*
 *  FUNCTION:  VECTOR_FLT_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data. 
 */
FLT VECTOR_FLT_Pop( VECTOR_FLT*   vec );

/*
 *  FUNCTION:  VECTOR_FLT_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
void VECTOR_FLT_Set(   VECTOR_FLT*     vec, 
                        int            idx, 
                        FLT            val );

/*
 *  FUNCTION:  VECTOR_FLT_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
FLT VECTOR_FLT_Get(  VECTOR_FLT*   vec, 
                     int           idx );

/*
 *  FUNCTION:  VECTOR_FLT_Get_Ref()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return pointer to data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
FLT* VECTOR_FLT_Get_Ref(   VECTOR_FLT*   vec, 
                           int           idx );

/*
 *  FUNCTION:  VECTOR_FLT_Search()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. Returns first instance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of first instance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_FLT_Search(  VECTOR_FLT*   vec, 
                        FLT           val );

/*
 *  FUNCTION:  VECTOR_FLT_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality, 
 *             pos if <vec_A> > <vec_B>,  
 *             neg if <vec_A> < <vec_B>.
 */
int VECTOR_FLT_Compare(    VECTOR_FLT*   vec_A, 
                           VECTOR_FLT*   vec_B );

/*
 *  FUNCTION:  VECTOR_FLT_Sort()
 *  SYNOPSIS:  Sort <vec> data array in ascending order.
 */
void VECTOR_FLT_Sort( VECTOR_FLT*    vec );

/*
 *  FUNCTION:  VECTOR_FLT_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
void VECTOR_FLT_Dump(  VECTOR_FLT*   vec,
                           FILE*             fp );

#endif 