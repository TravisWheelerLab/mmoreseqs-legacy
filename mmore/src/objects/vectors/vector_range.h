/*******************************************************************************
 *  FILE:      vector_range.c
 *  PURPOSE:   VECTOR_RANGE Object Functions
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _VECTOR_RANGE_H
#define _VECTOR_RANGE_H

/*! FUNCTION:  VECTOR_RANGE_Create()
 *  SYNOPSIS:  Allocates new VECTOR_RANGE object and returns pointer.
 */
VECTOR_RANGE* 
VECTOR_RANGE_Create();

/*! FUNCTION:  VECTOR_RANGE_Create_by_Size()
 *  SYNOPSIS:  Allocates new VECTOR_RANGE object at specific size and returns pointer.
 */
VECTOR_RANGE* 
VECTOR_RANGE_Create_by_Size( size_t    size );

/*! FUNCTION:  VECTOR_RANGE_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_RANGE and returns NULL pointer.
 */
VECTOR_RANGE* 
VECTOR_RANGE_Destroy( VECTOR_RANGE*   vec );

/*! FUNCTION:  VECTOR_RANGE_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_RANGE object by resetting size counter (no realloc) .
 */
STATUS_FLAG 
VECTOR_RANGE_Reuse( VECTOR_RANGE*   vec );

/*! FUNCTION:  VECTOR_RANGE_Fill()
 *  SYNOPSIS:  Fill VECTOR_RANGE object with val.
 */
STATUS_FLAG 
VECTOR_RANGE_Fill(  VECTOR_RANGE*   vec, 
                  RANGE           val );

/*! FUNCTION:  VECTOR_RANGE_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object. 
 *             Creates new VECTOR_RANGE for <dest> if <dest> is NULL.
 */
VECTOR_RANGE* 
VECTOR_RANGE_Copy(  VECTOR_RANGE*   dest, 
                  VECTOR_RANGE*   src );

/*! FUNCTION:  VECTOR_RANGE_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>. 
 */
STATUS_FLAG 
VECTOR_RANGE_Resize(   VECTOR_RANGE*    vec, 
                     size_t         size );

/*! FUNCTION:  VECTOR_RANGE_GetArray()
 *  SYNOPSIS:  Get <data> array from <vec>.
 */
RANGE* 
VECTOR_RANGE_GetArray(   VECTOR_RANGE*   vec );

/*! FUNCTION:  VECTOR_RANGE_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>. 
 */
STATUS_FLAG 
VECTOR_RANGE_GrowTo(   VECTOR_RANGE*    vec, 
                     size_t         size );

/*! FUNCTION:  VECTOR_RANGE_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
STATUS_FLAG 
VECTOR_RANGE_Push(  VECTOR_RANGE*   vec, 
                  RANGE           val );

/*! FUNCTION:  VECTOR_RANGE_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
STATUS_FLAG 
VECTOR_RANGE_Pushback(    VECTOR_RANGE*   vec, 
                        RANGE           val );

/*! FUNCTION:  VECTOR_RANGE_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data. 
 */
RANGE 
VECTOR_RANGE_Pop( VECTOR_RANGE*   vec );

/*! FUNCTION:  VECTOR_RANGE_Append()
 *  SYNOPSIS:  Push <append> data array of length <L> onto the end of <vec> data array. 
 */
STATUS_FLAG 
VECTOR_RANGE_Append(   VECTOR_RANGE*   vec, 
                     RANGE*          append,
                     size_t        L );

/*! FUNCTION:  VECTOR_RANGE_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
STATUS_FLAG 
VECTOR_RANGE_Set(   VECTOR_RANGE*       vec, 
                  int               idx, 
                  RANGE               val );

/*! FUNCTION:  VECTOR_RANGE_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
RANGE 
VECTOR_RANGE_Get(   VECTOR_RANGE*   vec, 
                  int           idx );

/*! FUNCTION:  VECTOR_RANGE_GetX()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return pointer to data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 *  RETURN:    Pointer to location to <vec> idx.
 */
RANGE* 
VECTOR_RANGE_GetX(    VECTOR_RANGE*   vec, 
                     int           idx );

/*! FUNCTION:  VECTOR_RANGE_GetSize()
 *  SYNOPSIS:  Get utilized length of <vec>.
 */
int 
VECTOR_RANGE_GetSize(   VECTOR_RANGE*   vec );

/*! FUNCTION:  VECTOR_RANGE_SetSize()
 *  SYNOPSIS:  Set utilized length of <vec>
 *  RETURN:    Pointer to location to <vec> idx.
 */
STATUS_FLAG 
VECTOR_RANGE_SetSize( VECTOR_RANGE*   vec, 
                     size_t        size );

/*! FUNCTION:  VECTOR_RANGE_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance found.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int 
VECTOR_RANGE_Search(   VECTOR_RANGE*   vec, 
                     RANGE           val );

/*! FUNCTION:  VECTOR_RANGE_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of first occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int 
VECTOR_RANGE_Search_First(   VECTOR_RANGE*   vec, 
                           RANGE           val );

/*! FUNCTION:  VECTOR_RANGE_Search_Last()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the last occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of last occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int 
VECTOR_RANGE_Search_Last( VECTOR_RANGE*   vec, 
                        RANGE           val );

/*! FUNCTION:  VECTOR_RANGE_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality, 
 *             pos if <vec_A> > <vec_B>,  
 *             neg if <vec_A> < <vec_B>.
 */
int 
VECTOR_RANGE_Compare(  VECTOR_RANGE*   vec_A, 
                     VECTOR_RANGE*   vec_B );

/*! FUNCTION:  VECTOR_RANGE_Sort()
 *  SYNOPSIS:  Sorts <vec> data array in ascending order. In-place.
 */
STATUS_FLAG 
VECTOR_RANGE_Sort( VECTOR_RANGE*    vec );

/*! FUNCTION:  VECTOR_RANGE_Sort_Sub()
 *  SYNOPSIS:  Subcall - Sorts <vec> data array in ascending order on range (beg,end]. 
 *             Uses quicksort until length of subarray falls below threshold, then selection sort.
 */
STATUS_FLAG 
VECTOR_RANGE_Sort_Sub(    VECTOR_RANGE*    vec,
                        int            beg,
                        int            end );

/*! FUNCTION:  VECTOR_RANGE_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs selection sort on <vec> data array in ascending order on range (beg,end]. 
 */
STATUS_FLAG 
VECTOR_RANGE_Sort_Sub_Selectsort(  VECTOR_RANGE*    vec,
                                 int            beg,
                                 int            end );

/*! FUNCTION:  VECTOR_RANGE_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs quicksort on <vec> data array in ascending order on range (beg,end]. 
 */
STATUS_FLAG 
VECTOR_RANGE_Sort_Sub_Quicksort(   VECTOR_RANGE*    vec,
                                 int            beg,
                                 int            end );

/*! FUNCTION:  VECTOR_RANGE_Swap()
 *  SYNOPSIS:  Swaps the values of <vec> at indexes <i> and <j>.
 *             Warning: Only checks for Out-of-Bounds when in DEBUG.
 */
STATUS_FLAG 
VECTOR_RANGE_Swap(  VECTOR_RANGE*    vec,
                  int            i,
                  int            j );

/*! FUNCTION:  VECTOR_RANGE_Reverse()
 *  SYNOPSIS:  Reverse the ordering of array.
 */
STATUS_FLAG 
VECTOR_RANGE_Reverse(   VECTOR_RANGE*    vec );

/*! FUNCTION:  VECTOR_RANGE_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
STATUS_FLAG 
VECTOR_RANGE_Dump(  VECTOR_RANGE*    vec,
                  FILE*          fp );

/*! FUNCTION:  VECTOR_RANGE_Unit_Test()
 *  SYNOPSIS:  Perform unit test for VECTOR_RANGE.
 */
STATUS_FLAG 
VECTOR_RANGE_Unit_Test();

#endif 