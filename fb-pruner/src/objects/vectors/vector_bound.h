/*******************************************************************************
 *  FILE:      vector_bound.c
 *  PURPOSE:   VECTOR_BOUND Object Functions
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _VECTOR_BOUND_H
#define _VECTOR_BOUND_H

/*
 *  FUNCTION:  VECTOR_BOUND_Create()
 *  SYNOPSIS:  Allocates new VECTOR_BOUND object and returns pointer.
 */
VECTOR_BOUND* VECTOR_BOUND_Create();

/*
 *  FUNCTION:  VECTOR_BOUND_Create()
 *  SYNOPSIS:  Allocates new VECTOR_BOUND object at specific size and returns pointer.
 */
VECTOR_BOUND* VECTOR_BOUND_Create_by_Size( int    size );

/*
 *  FUNCTION:  VECTOR_BOUND_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_BOUND.
 */
void* VECTOR_BOUND_Destroy( VECTOR_BOUND*   vec );

/*
 *  FUNCTION:  VECTOR_BOUND_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_BOUND object by resetting size counter (no realloc) .
 */
void VECTOR_BOUND_Reuse( VECTOR_BOUND*   vec );

/*
 *  FUNCTION:  VECTOR_BOUND_Fill()
 *  SYNOPSIS:  Fill VECTOR_BOUND object with val.
 */
void VECTOR_BOUND_Fill(   VECTOR_BOUND*   vec, 
                        BOUND           val );

/*
 *  FUNCTION:  VECTOR_BOUND_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object. 
 *             Creates new VECTOR_BOUND for <dest> if <dest> is NULL.
 */
VECTOR_BOUND* VECTOR_BOUND_Copy(  VECTOR_BOUND*   src, 
                              VECTOR_BOUND*   dest );

/*
 *  FUNCTION:  VECTOR_BOUND_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>. 
 */
void VECTOR_BOUND_Resize(    VECTOR_BOUND*    vec, 
                           int            size );

/*
 *  FUNCTION:  VECTOR_BOUND_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>. 
 */
void VECTOR_BOUND_GrowTo(   VECTOR_BOUND*     vec, 
                           int            size );

/*
 *  FUNCTION:  VECTOR_BOUND_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
void VECTOR_BOUND_Push(   VECTOR_BOUND*   vec, 
                        BOUND           val );

/*
 *  FUNCTION:  VECTOR_BOUND_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
void VECTOR_BOUND_Pushback(  VECTOR_BOUND*   vec, 
                           BOUND           val );

/*
 *  FUNCTION:  VECTOR_BOUND_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data. 
 */
BOUND VECTOR_BOUND_Pop( VECTOR_BOUND*   vec );

/*
 *  FUNCTION:  VECTOR_BOUND_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
void VECTOR_BOUND_Set(   VECTOR_BOUND*     vec, 
                        int            idx, 
                        BOUND            val );

/*
 *  FUNCTION:  VECTOR_BOUND_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
BOUND VECTOR_BOUND_Get(    VECTOR_BOUND*   vec, 
                        int           idx );

/*
 *  FUNCTION:  VECTOR_BOUND_Get_X()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return pointer to data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 *  RETURN:    Pointer to location to <vec> idx.
 */
BOUND* VECTOR_BOUND_Get_X(  VECTOR_BOUND*   vec, 
                        int           idx );

/*
 *  FUNCTION:  VECTOR_BOUND_Get_Size()
 *  SYNOPSIS:  Get utilized length of <vec>.
 */
int VECTOR_BOUND_Get_Size(   VECTOR_BOUND*   vec );

/*
 *  FUNCTION:  VECTOR_BOUND_Set_Size()
 *  SYNOPSIS:  Set utilized length of <vec>
 *  RETURN:    Pointer to location to <vec> idx.
 */
void VECTOR_BOUND_Set_Size(  VECTOR_BOUND*   vec, 
                           int           size );

/*
 *  FUNCTION:  VECTOR_BOUND_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance found.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_BOUND_Search(  VECTOR_BOUND*   vec, 
                        BOUND           val );

/*
 *  FUNCTION:  VECTOR_BOUND_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of first occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_BOUND_Search_First(  VECTOR_BOUND*   vec, 
                              BOUND           val );

/*
 *  FUNCTION:  VECTOR_BOUND_Search_Last()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the last occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of last occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_BOUND_Search_Last(   VECTOR_BOUND*   vec, 
                              BOUND           val );

/*
 *  FUNCTION:  VECTOR_BOUND_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality, 
 *             pos if <vec_A> > <vec_B>,  
 *             neg if <vec_A> < <vec_B>.
 */
int VECTOR_BOUND_Compare(    VECTOR_BOUND*   vec_A, 
                           VECTOR_BOUND*   vec_B );

/*
 *  FUNCTION:  VECTOR_BOUND_Sort()
 *  SYNOPSIS:  Sorts <vec> data array in ascending order. In-place.
 */
void VECTOR_BOUND_Sort( VECTOR_BOUND*    vec );

/*
 *  FUNCTION:  VECTOR_BOUND_Sort_Sub()
 *  SYNOPSIS:  Subcall - Sorts <vec> data array in ascending order on range (beg,end]. 
 *             Uses quicksort until length of subarray falls below threshold, then selection sort.
 */
void VECTOR_BOUND_Sort_Sub(  VECTOR_BOUND*    vec,
                           int            beg,
                           int            end );

/*
 *  FUNCTION:  VECTOR_BOUND_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs selection sort on <vec> data array in ascending order on range (beg,end]. 
 */
void VECTOR_BOUND_Sort_Sub_Selectsort(   VECTOR_BOUND*    vec,
                                       int            beg,
                                       int            end );

/*
 *  FUNCTION:  VECTOR_BOUND_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs quicksort on <vec> data array in ascending order on range (beg,end]. 
 */
void VECTOR_BOUND_Sort_Sub_Quicksort( VECTOR_BOUND*    vec,
                                    int            beg,
                                    int            end );

/*
 *  FUNCTION:  VECTOR_BOUND_Swap()
 *  SYNOPSIS:  Swaps the values of <vec> at indexes <i> and <j>.
 *             Warning: Only checks for Out-of-Bounds when in DEBUG.
 */
void VECTOR_BOUND_Swap(   VECTOR_BOUND*    vec,
                        int            i,
                        int            j );

/*
 *  FUNCTION:  VECTOR_BOUND_Reverse()
 *  SYNOPSIS:  Reverse the ordering of array.
 */
void VECTOR_BOUND_Reverse(   VECTOR_BOUND*    vec );

/*
 *  FUNCTION:  VECTOR_BOUND_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
void VECTOR_BOUND_Dump(   VECTOR_BOUND*    vec,
                        FILE*          fp );

/*
 *  FUNCTION:  VECTOR_BOUND_Unit_Test()
 *  SYNOPSIS:  Perform unit test for VECTOR_BOUND.
 */
void VECTOR_BOUND_Unit_Test();

#endif 