/*******************************************************************************
 *  FILE:      vector_int.c
 *  PURPOSE:   VECTOR_INT Object Functions
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _VECTOR_INT_H
#define _VECTOR_INT_H

/*
 *  FUNCTION:  VECTOR_INT_Create()
 *  SYNOPSIS:  Allocates new VECTOR_INT object and returns pointer.
 */
VECTOR_INT* VECTOR_INT_Create();

/*
 *  FUNCTION:  VECTOR_INT_Create()
 *  SYNOPSIS:  Allocates new VECTOR_INT object at specific size and returns pointer.
 */
VECTOR_INT* VECTOR_INT_Create_by_Size( int    size );

/*
 *  FUNCTION:  VECTOR_INT_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_INT.
 */
void* VECTOR_INT_Destroy( VECTOR_INT*   vec );

/*
 *  FUNCTION:  VECTOR_INT_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_INT object by resetting size counter (no realloc) .
 */
void VECTOR_INT_Reuse( VECTOR_INT*   vec );

/*
 *  FUNCTION:  VECTOR_INT_Fill()
 *  SYNOPSIS:  Fill VECTOR_INT object with val.
 */
void VECTOR_INT_Fill(   VECTOR_INT*   vec, 
                        INT           val );

/*
 *  FUNCTION:  VECTOR_INT_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object. 
 *             Creates new VECTOR_INT for <dest> if <dest> is NULL.
 */
VECTOR_INT* VECTOR_INT_Copy(  VECTOR_INT*   src, 
                              VECTOR_INT*   dest );

/*
 *  FUNCTION:  VECTOR_INT_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>. 
 */
void VECTOR_INT_Resize(    VECTOR_INT*    vec, 
                           int            size );

/*
 *  FUNCTION:  VECTOR_INT_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>. 
 */
void VECTOR_INT_GrowTo(   VECTOR_INT*     vec, 
                           int            size );

/*
 *  FUNCTION:  VECTOR_INT_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
void VECTOR_INT_Push(   VECTOR_INT*   vec, 
                        INT           val );

/*
 *  FUNCTION:  VECTOR_INT_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
void VECTOR_INT_Pushback(  VECTOR_INT*   vec, 
                           INT           val );

/*
 *  FUNCTION:  VECTOR_INT_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data. 
 */
INT VECTOR_INT_Pop( VECTOR_INT*   vec );

/*
 *  FUNCTION:  VECTOR_INT_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
void VECTOR_INT_Set(   VECTOR_INT*     vec, 
                        int            idx, 
                        INT            val );

/*
 *  FUNCTION:  VECTOR_INT_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
INT VECTOR_INT_Get(    VECTOR_INT*   vec, 
                        int           idx );

/*
 *  FUNCTION:  VECTOR_INT_Get_X()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return pointer to data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 *  RETURN:    Pointer to location to <vec> idx.
 */
INT* VECTOR_INT_Get_X(  VECTOR_INT*   vec, 
                        int           idx );

/*
 *  FUNCTION:  VECTOR_INT_Get_Size()
 *  SYNOPSIS:  Get utilized length of <vec>.
 */
INT VECTOR_INT_Get_Size(   VECTOR_INT*   vec );

/*
 *  FUNCTION:  VECTOR_INT_Set_Size()
 *  SYNOPSIS:  Set utilized length of <vec>
 *  RETURN:    Pointer to location to <vec> idx.
 */
void VECTOR_INT_Set_Size(  VECTOR_INT*   vec, 
                           int           size );

/*
 *  FUNCTION:  VECTOR_INT_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance found.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_INT_Search(  VECTOR_INT*   vec, 
                        INT           val );

/*
 *  FUNCTION:  VECTOR_INT_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of first occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_INT_Search_First(  VECTOR_INT*   vec, 
                              INT           val );

/*
 *  FUNCTION:  VECTOR_INT_Search_Last()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the last occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of last occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_INT_Search_Last(   VECTOR_INT*   vec, 
                              INT           val );

/*
 *  FUNCTION:  VECTOR_INT_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality, 
 *             pos if <vec_A> > <vec_B>,  
 *             neg if <vec_A> < <vec_B>.
 */
int VECTOR_INT_Compare(    VECTOR_INT*   vec_A, 
                           VECTOR_INT*   vec_B );

/*
 *  FUNCTION:  VECTOR_INT_Sort()
 *  SYNOPSIS:  Sorts <vec> data array in ascending order. In-place.
 */
void VECTOR_INT_Sort( VECTOR_INT*    vec );

/*
 *  FUNCTION:  VECTOR_INT_Sort_Sub()
 *  SYNOPSIS:  Subcall - Sorts <vec> data array in ascending order on range (beg,end]. 
 *             Uses quicksort until length of subarray falls below threshold, then selection sort.
 */
void VECTOR_INT_Sort_Sub(  VECTOR_INT*    vec,
                           int            beg,
                           int            end );

/*
 *  FUNCTION:  VECTOR_INT_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs selection sort on <vec> data array in ascending order on range (beg,end]. 
 */
void VECTOR_INT_Sort_Sub_Selectsort(   VECTOR_INT*    vec,
                                       int            beg,
                                       int            end );

/*
 *  FUNCTION:  VECTOR_INT_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs quicksort on <vec> data array in ascending order on range (beg,end]. 
 */
void VECTOR_INT_Sort_Sub_Quicksort( VECTOR_INT*    vec,
                                    int            beg,
                                    int            end );

/*
 *  FUNCTION:  VECTOR_INT_Swap()
 *  SYNOPSIS:  Swaps the values of <vec> at indexes <i> and <j>.
 *             Warning: Only checks for Out-of-Bounds when in DEBUG.
 */
void VECTOR_INT_Swap(   VECTOR_INT*    vec,
                        int            i,
                        int            j );

/*
 *  FUNCTION:  VECTOR_INT_Reverse()
 *  SYNOPSIS:  Reverse the ordering of array.
 */
void VECTOR_INT_Reverse(   VECTOR_INT*    vec );

/*
 *  FUNCTION:  VECTOR_INT_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
void VECTOR_INT_Dump(   VECTOR_INT*    vec,
                        FILE*          fp );

#endif 