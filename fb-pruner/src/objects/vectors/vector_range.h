/*******************************************************************************
 *  FILE:      vector_range.c
 *  PURPOSE:   VECTOR_RANGE Object Functions
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _VECTOR_RANGE_H
#define _VECTOR_RANGE_H

/*
 *  FUNCTION:  VECTOR_RANGE_Create()
 *  SYNOPSIS:  Allocates new VECTOR_RANGE object and returns pointer.
 */
VECTOR_RANGE* VECTOR_RANGE_Create();

/*
 *  FUNCTION:  VECTOR_RANGE_Create()
 *  SYNOPSIS:  Allocates new VECTOR_RANGE object at specific size and returns pointer.
 */
VECTOR_RANGE* VECTOR_RANGE_Create_by_Size( int    size );

/*
 *  FUNCTION:  VECTOR_RANGE_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_RANGE.
 */
void* VECTOR_RANGE_Destroy( VECTOR_RANGE*   vec );

/*
 *  FUNCTION:  VECTOR_RANGE_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_RANGE object by resetting size counter (no realloc) .
 */
void VECTOR_RANGE_Reuse( VECTOR_RANGE*   vec );

/*
 *  FUNCTION:  VECTOR_RANGE_Fill()
 *  SYNOPSIS:  Fill VECTOR_RANGE object with val.
 */
void VECTOR_RANGE_Fill(   VECTOR_RANGE*   vec, 
                        RANGE           val );

/*
 *  FUNCTION:  VECTOR_RANGE_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object. 
 *             Creates new VECTOR_RANGE for <dest> if <dest> is NULL.
 */
VECTOR_RANGE* VECTOR_RANGE_Copy(  VECTOR_RANGE*   src, 
                              VECTOR_RANGE*   dest );

/*
 *  FUNCTION:  VECTOR_RANGE_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>. 
 */
void VECTOR_RANGE_Resize(    VECTOR_RANGE*    vec, 
                           int            size );

/*
 *  FUNCTION:  VECTOR_RANGE_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>. 
 */
void VECTOR_RANGE_GrowTo(   VECTOR_RANGE*     vec, 
                           int            size );

/*
 *  FUNCTION:  VECTOR_RANGE_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
void VECTOR_RANGE_Push(   VECTOR_RANGE*   vec, 
                        RANGE           val );

/*
 *  FUNCTION:  VECTOR_RANGE_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
void VECTOR_RANGE_Pushback(  VECTOR_RANGE*   vec, 
                           RANGE           val );

/*
 *  FUNCTION:  VECTOR_RANGE_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data. 
 */
RANGE VECTOR_RANGE_Pop( VECTOR_RANGE*   vec );

/*
 *  FUNCTION:  VECTOR_RANGE_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
void VECTOR_RANGE_Set(   VECTOR_RANGE*     vec, 
                        int            idx, 
                        RANGE            val );

/*
 *  FUNCTION:  VECTOR_RANGE_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
RANGE VECTOR_RANGE_Get(    VECTOR_RANGE*   vec, 
                        int           idx );

/*
 *  FUNCTION:  VECTOR_RANGE_Get_X()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return pointer to data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 *  RETURN:    Pointer to location to <vec> idx.
 */
RANGE* VECTOR_RANGE_Get_X(  VECTOR_RANGE*   vec, 
                        int           idx );

/*
 *  FUNCTION:  VECTOR_RANGE_Get_Size()
 *  SYNOPSIS:  Get utilized length of <vec>.
 */
RANGE VECTOR_RANGE_Get_Size(   VECTOR_RANGE*   vec );

/*
 *  FUNCTION:  VECTOR_RANGE_Set_Size()
 *  SYNOPSIS:  Set utilized length of <vec>
 *  RETURN:    Pointer to location to <vec> idx.
 */
void VECTOR_RANGE_Set_Size(  VECTOR_RANGE*   vec, 
                           int           size );

/*
 *  FUNCTION:  VECTOR_RANGE_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance found.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_RANGE_Search(  VECTOR_RANGE*   vec, 
                        RANGE           val );

/*
 *  FUNCTION:  VECTOR_RANGE_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of first occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_RANGE_Search_First(  VECTOR_RANGE*   vec, 
                              RANGE           val );

/*
 *  FUNCTION:  VECTOR_RANGE_Search_Last()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the last occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of last occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_RANGE_Search_Last(   VECTOR_RANGE*   vec, 
                              RANGE           val );

/*
 *  FUNCTION:  VECTOR_RANGE_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality, 
 *             pos if <vec_A> > <vec_B>,  
 *             neg if <vec_A> < <vec_B>.
 */
int VECTOR_RANGE_Compare(    VECTOR_RANGE*   vec_A, 
                           VECTOR_RANGE*   vec_B );

/*
 *  FUNCTION:  VECTOR_RANGE_Sort()
 *  SYNOPSIS:  Sorts <vec> data array in ascending order. In-place.
 */
void VECTOR_RANGE_Sort( VECTOR_RANGE*    vec );

/*
 *  FUNCTION:  VECTOR_RANGE_Sort_Sub()
 *  SYNOPSIS:  Subcall - Sorts <vec> data array in ascending order on range (beg,end]. 
 *             Uses quicksort until length of subarray falls below threshold, then selection sort.
 */
void VECTOR_RANGE_Sort_Sub(  VECTOR_RANGE*    vec,
                           int            beg,
                           int            end );

/*
 *  FUNCTION:  VECTOR_RANGE_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs selection sort on <vec> data array in ascending order on range (beg,end]. 
 */
void VECTOR_RANGE_Sort_Sub_Selectsort(   VECTOR_RANGE*    vec,
                                       int            beg,
                                       int            end );

/*
 *  FUNCTION:  VECTOR_RANGE_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs quicksort on <vec> data array in ascending order on range (beg,end]. 
 */
void VECTOR_RANGE_Sort_Sub_Quicksort( VECTOR_RANGE*    vec,
                                    int            beg,
                                    int            end );

/*
 *  FUNCTION:  VECTOR_RANGE_Swap()
 *  SYNOPSIS:  Swaps the values of <vec> at indexes <i> and <j>.
 *             Warning: Only checks for Out-of-Bounds when in DEBUG.
 */
void VECTOR_RANGE_Swap(   VECTOR_RANGE*    vec,
                        int            i,
                        int            j );

/*
 *  FUNCTION:  VECTOR_RANGE_Reverse()
 *  SYNOPSIS:  Reverse the ordering of array.
 */
void VECTOR_RANGE_Reverse(   VECTOR_RANGE*    vec );

/*
 *  FUNCTION:  VECTOR_RANGE_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
void VECTOR_RANGE_Dump(   VECTOR_RANGE*    vec,
                        FILE*          fp );

#endif 