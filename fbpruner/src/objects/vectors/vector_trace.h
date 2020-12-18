/*******************************************************************************
 *  FILE:      vector_trace.c
 *  PURPOSE:   VECTOR_TRACE Object Functions
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _VECTOR_TRACE_H
#define _VECTOR_TRACE_H

/*
 *  FUNCTION:  VECTOR_TRACE_Create()
 *  SYNOPSIS:  Allocates new VECTOR_TRACE object and returns pointer.
 */
VECTOR_TRACE* VECTOR_TRACE_Create();

/*
 *  FUNCTION:  VECTOR_TRACE_Create()
 *  SYNOPSIS:  Allocates new VECTOR_TRACE object at specific size and returns pointer.
 */
VECTOR_TRACE* VECTOR_TRACE_Create_by_Size( int    size );

/*
 *  FUNCTION:  VECTOR_TRACE_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_TRACE.
 */
void* VECTOR_TRACE_Destroy( VECTOR_TRACE*   vec );

/*
 *  FUNCTION:  VECTOR_TRACE_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_TRACE object by resetting size counter (no realloc) .
 */
void VECTOR_TRACE_Reuse( VECTOR_TRACE*   vec );

/*
 *  FUNCTION:  VECTOR_TRACE_Fill()
 *  SYNOPSIS:  Fill VECTOR_TRACE object with val.
 */
void VECTOR_TRACE_Fill(   VECTOR_TRACE*   vec, 
                        TRACE           val );

/*
 *  FUNCTION:  VECTOR_TRACE_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object. 
 *             Creates new VECTOR_TRACE for <dest> if <dest> is NULL.
 */
VECTOR_TRACE* VECTOR_TRACE_Copy(  VECTOR_TRACE*   src, 
                              VECTOR_TRACE*   dest );

/*
 *  FUNCTION:  VECTOR_TRACE_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>. 
 */
void VECTOR_TRACE_Resize(    VECTOR_TRACE*    vec, 
                           int            size );

/*
 *  FUNCTION:  VECTOR_TRACE_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>. 
 */
void VECTOR_TRACE_GrowTo(   VECTOR_TRACE*     vec, 
                           int            size );

/*
 *  FUNCTION:  VECTOR_TRACE_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
void VECTOR_TRACE_Push(   VECTOR_TRACE*   vec, 
                        TRACE           val );

/*
 *  FUNCTION:  VECTOR_TRACE_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
void VECTOR_TRACE_Pushback(  VECTOR_TRACE*   vec, 
                           TRACE           val );

/*
 *  FUNCTION:  VECTOR_TRACE_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data. 
 */
TRACE VECTOR_TRACE_Pop( VECTOR_TRACE*   vec );

/*
 *  FUNCTION:  VECTOR_TRACE_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
void VECTOR_TRACE_Set(   VECTOR_TRACE*     vec, 
                        int            idx, 
                        TRACE            val );

/*
 *  FUNCTION:  VECTOR_TRACE_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
TRACE VECTOR_TRACE_Get(    VECTOR_TRACE*   vec, 
                        int           idx );

/*
 *  FUNCTION:  VECTOR_TRACE_Get_X()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return pointer to data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 *  RETURN:    Pointer to location to <vec> idx.
 */
TRACE* VECTOR_TRACE_Get_X(  VECTOR_TRACE*   vec, 
                        int           idx );

/*
 *  FUNCTION:  VECTOR_TRACE_Get_Size()
 *  SYNOPSIS:  Get utilized length of <vec>.
 */
int VECTOR_TRACE_Get_Size(   VECTOR_TRACE*   vec );

/*
 *  FUNCTION:  VECTOR_TRACE_Set_Size()
 *  SYNOPSIS:  Set utilized length of <vec>
 *  RETURN:    Pointer to location to <vec> idx.
 */
void VECTOR_TRACE_Set_Size(  VECTOR_TRACE*   vec, 
                           int           size );

/*
 *  FUNCTION:  VECTOR_TRACE_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance found.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_TRACE_Search(  VECTOR_TRACE*   vec, 
                        TRACE           val );

/*
 *  FUNCTION:  VECTOR_TRACE_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of first occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_TRACE_Search_First(  VECTOR_TRACE*   vec, 
                              TRACE           val );

/*
 *  FUNCTION:  VECTOR_TRACE_Search_Last()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the last occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of last occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_TRACE_Search_Last(   VECTOR_TRACE*   vec, 
                              TRACE           val );

/*
 *  FUNCTION:  VECTOR_TRACE_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality, 
 *             pos if <vec_A> > <vec_B>,  
 *             neg if <vec_A> < <vec_B>.
 */
int VECTOR_TRACE_Compare(    VECTOR_TRACE*   vec_A, 
                           VECTOR_TRACE*   vec_B );

/*
 *  FUNCTION:  VECTOR_TRACE_Sort()
 *  SYNOPSIS:  Sorts <vec> data array in ascending order. In-place.
 */
void VECTOR_TRACE_Sort( VECTOR_TRACE*    vec );

/*
 *  FUNCTION:  VECTOR_TRACE_Sort_Sub()
 *  SYNOPSIS:  Subcall - Sorts <vec> data array in ascending order on range (beg,end]. 
 *             Uses quicksort until length of subarray falls below threshold, then selection sort.
 */
void VECTOR_TRACE_Sort_Sub(  VECTOR_TRACE*    vec,
                           int            beg,
                           int            end );

/*
 *  FUNCTION:  VECTOR_TRACE_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs selection sort on <vec> data array in ascending order on range (beg,end]. 
 */
void VECTOR_TRACE_Sort_Sub_Selectsort(   VECTOR_TRACE*    vec,
                                       int            beg,
                                       int            end );

/*
 *  FUNCTION:  VECTOR_TRACE_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs quicksort on <vec> data array in ascending order on range (beg,end]. 
 */
void VECTOR_TRACE_Sort_Sub_Quicksort( VECTOR_TRACE*    vec,
                                    int            beg,
                                    int            end );

/*
 *  FUNCTION:  VECTOR_TRACE_Swap()
 *  SYNOPSIS:  Swaps the values of <vec> at indexes <i> and <j>.
 *             Warning: Only checks for Out-of-Bounds when in DEBUG.
 */
void VECTOR_TRACE_Swap(   VECTOR_TRACE*    vec,
                        int            i,
                        int            j );

/*
 *  FUNCTION:  VECTOR_TRACE_Reverse()
 *  SYNOPSIS:  Reverse the ordering of array.
 */
void VECTOR_TRACE_Reverse(   VECTOR_TRACE*    vec );

/*
 *  FUNCTION:  VECTOR_TRACE_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
void VECTOR_TRACE_Dump(   VECTOR_TRACE*    vec,
                        FILE*          fp );

/*
 *  FUNCTION:  VECTOR_TRACE_Unit_Test()
 *  SYNOPSIS:  Perform unit test for VECTOR_TRACE.
 */
void VECTOR_TRACE_Unit_Test();

#endif 