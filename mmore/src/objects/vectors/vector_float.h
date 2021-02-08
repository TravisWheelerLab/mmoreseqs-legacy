/*******************************************************************************
 *  FILE:      vector_float.c
 *  PURPOSE:   VECTOR_FLT Object Functions
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _VECTOR_FLT_H
#define _VECTOR_FLT_H

/*! FUNCTION:  VECTOR_FLT_Create()
 *  SYNOPSIS:  Allocates new VECTOR_FLT object and returns pointer.
 */
VECTOR_FLT* 
VECTOR_FLT_Create();

/*! FUNCTION:  VECTOR_FLT_Create_by_Size()
 *  SYNOPSIS:  Allocates new VECTOR_FLT object at specific size and returns pointer.
 */
VECTOR_FLT* 
VECTOR_FLT_Create_by_Size( size_t    size );

/*! FUNCTION:  VECTOR_FLT_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_FLT and returns NULL pointer.
 */
VECTOR_FLT* 
VECTOR_FLT_Destroy( VECTOR_FLT*   vec );

/*! FUNCTION:  VECTOR_FLT_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_FLT object by resetting size counter (no realloc) .
 */
STATUS_FLAG 
VECTOR_FLT_Reuse( VECTOR_FLT*   vec );

/*! FUNCTION:  VECTOR_FLT_Fill()
 *  SYNOPSIS:  Fill VECTOR_FLT object with val.
 */
STATUS_FLAG 
VECTOR_FLT_Fill(  VECTOR_FLT*   vec, 
                  FLT           val );

/*! FUNCTION:  VECTOR_FLT_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object. 
 *             Creates new VECTOR_FLT for <dest> if <dest> is NULL.
 */
VECTOR_FLT* 
VECTOR_FLT_Copy(  VECTOR_FLT*   dest, 
                  VECTOR_FLT*   src );

/*! FUNCTION:  VECTOR_FLT_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>. 
 */
STATUS_FLAG 
VECTOR_FLT_Resize(   VECTOR_FLT*    vec, 
                     size_t         size );

/*! FUNCTION:  VECTOR_FLT_GetArray()
 *  SYNOPSIS:  Get <data> array from <vec>.
 */
FLT* 
VECTOR_FLT_GetArray(   VECTOR_FLT*   vec );

/*! FUNCTION:  VECTOR_FLT_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>. 
 */
STATUS_FLAG 
VECTOR_FLT_GrowTo(   VECTOR_FLT*    vec, 
                     size_t         size );

/*! FUNCTION:  VECTOR_FLT_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
STATUS_FLAG 
VECTOR_FLT_Push(  VECTOR_FLT*   vec, 
                  FLT           val );

/*! FUNCTION:  VECTOR_FLT_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
STATUS_FLAG 
VECTOR_FLT_Pushback(    VECTOR_FLT*   vec, 
                        FLT           val );

/*! FUNCTION:  VECTOR_FLT_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data. 
 */
FLT 
VECTOR_FLT_Pop( VECTOR_FLT*   vec );

/*! FUNCTION:  VECTOR_FLT_Append()
 *  SYNOPSIS:  Push <append> data array of length <L> onto the end of <vec> data array. 
 */
STATUS_FLAG 
VECTOR_FLT_Append(   VECTOR_FLT*   vec, 
                     FLT*          append,
                     size_t        L );

/*! FUNCTION:  VECTOR_FLT_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
STATUS_FLAG 
VECTOR_FLT_Set(   VECTOR_FLT*       vec, 
                  int               idx, 
                  FLT               val );

/*! FUNCTION:  VECTOR_FLT_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
FLT 
VECTOR_FLT_Get(   VECTOR_FLT*   vec, 
                  int           idx );

/*! FUNCTION:  VECTOR_FLT_GetX()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return pointer to data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 *  RETURN:    Pointer to location to <vec> idx.
 */
FLT* 
VECTOR_FLT_GetX(    VECTOR_FLT*   vec, 
                     int           idx );

/*! FUNCTION:  VECTOR_FLT_GetSize()
 *  SYNOPSIS:  Get utilized length of <vec>.
 */
int 
VECTOR_FLT_GetSize(   VECTOR_FLT*   vec );

/*! FUNCTION:  VECTOR_FLT_SetSize()
 *  SYNOPSIS:  Set utilized length of <vec>
 *  RETURN:    Pointer to location to <vec> idx.
 */
STATUS_FLAG 
VECTOR_FLT_SetSize( VECTOR_FLT*   vec, 
                     size_t        size );

/*! FUNCTION:  VECTOR_FLT_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance found.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int 
VECTOR_FLT_Search(   VECTOR_FLT*   vec, 
                     FLT           val );

/*! FUNCTION:  VECTOR_FLT_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of first occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int 
VECTOR_FLT_Search_First(   VECTOR_FLT*   vec, 
                           FLT           val );

/*! FUNCTION:  VECTOR_FLT_Search_Last()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the last occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of last occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int 
VECTOR_FLT_Search_Last( VECTOR_FLT*   vec, 
                        FLT           val );

/*! FUNCTION:  VECTOR_FLT_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality, 
 *             pos if <vec_A> > <vec_B>,  
 *             neg if <vec_A> < <vec_B>.
 */
int 
VECTOR_FLT_Compare(  VECTOR_FLT*   vec_A, 
                     VECTOR_FLT*   vec_B );

/*! FUNCTION:  VECTOR_FLT_Sort()
 *  SYNOPSIS:  Sorts <vec> data array in ascending order. In-place.
 */
STATUS_FLAG 
VECTOR_FLT_Sort( VECTOR_FLT*    vec );

/*! FUNCTION:  VECTOR_FLT_Sort_Sub()
 *  SYNOPSIS:  Subcall - Sorts <vec> data array in ascending order on range (beg,end]. 
 *             Uses quicksort until length of subarray falls below threshold, then selection sort.
 */
STATUS_FLAG 
VECTOR_FLT_Sort_Sub(    VECTOR_FLT*    vec,
                        int            beg,
                        int            end );

/*! FUNCTION:  VECTOR_FLT_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs selection sort on <vec> data array in ascending order on range (beg,end]. 
 */
STATUS_FLAG 
VECTOR_FLT_Sort_Sub_Selectsort(  VECTOR_FLT*    vec,
                                 int            beg,
                                 int            end );

/*! FUNCTION:  VECTOR_FLT_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs quicksort on <vec> data array in ascending order on range (beg,end]. 
 */
STATUS_FLAG 
VECTOR_FLT_Sort_Sub_Quicksort(   VECTOR_FLT*    vec,
                                 int            beg,
                                 int            end );

/*! FUNCTION:  VECTOR_FLT_Swap()
 *  SYNOPSIS:  Swaps the values of <vec> at indexes <i> and <j>.
 *             Warning: Only checks for Out-of-Bounds when in DEBUG.
 */
STATUS_FLAG 
VECTOR_FLT_Swap(  VECTOR_FLT*    vec,
                  int            i,
                  int            j );

/*! FUNCTION:  VECTOR_FLT_Reverse()
 *  SYNOPSIS:  Reverse the ordering of array.
 */
STATUS_FLAG 
VECTOR_FLT_Reverse(   VECTOR_FLT*    vec );

/*! FUNCTION:  VECTOR_FLT_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
STATUS_FLAG 
VECTOR_FLT_Dump(  VECTOR_FLT*    vec,
                  FILE*          fp );

/*! FUNCTION:  VECTOR_FLT_Unit_Test()
 *  SYNOPSIS:  Perform unit test for VECTOR_FLT.
 */
STATUS_FLAG 
VECTOR_FLT_Unit_Test();

#endif 