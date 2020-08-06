/*******************************************************************************
 *  FILE:      vector_double.c
 *  PURPOSE:   VECTOR_DBL Object Functions
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _VECTOR_DBL_H
#define _VECTOR_DBL_H

/*
 *  FUNCTION:  VECTOR_DBL_Create()
 *  SYNOPSIS:  Allocates new VECTOR_DBL object and returns pointer.
 */
VECTOR_DBL* VECTOR_DBL_Create();

/*
 *  FUNCTION:  VECTOR_DBL_Create()
 *  SYNOPSIS:  Allocates new VECTOR_DBL object at specific size and returns pointer.
 */
VECTOR_DBL* VECTOR_DBL_Create_by_Size( int    size );

/*
 *  FUNCTION:  VECTOR_DBL_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_DBL.
 */
void* VECTOR_DBL_Destroy( VECTOR_DBL*   vec );

/*
 *  FUNCTION:  VECTOR_DBL_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_DBL object by resetting size counter (no realloc) .
 */
void VECTOR_DBL_Reuse( VECTOR_DBL*   vec );

/*
 *  FUNCTION:  VECTOR_DBL_Fill()
 *  SYNOPSIS:  Fill VECTOR_DBL object with val.
 */
void VECTOR_DBL_Fill(   VECTOR_DBL*   vec, 
                        DBL           val );

/*
 *  FUNCTION:  VECTOR_DBL_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object. 
 *             Creates new VECTOR_DBL for <dest> if <dest> is NULL.
 */
VECTOR_DBL* VECTOR_DBL_Copy(  VECTOR_DBL*   src, 
                              VECTOR_DBL*   dest );

/*
 *  FUNCTION:  VECTOR_DBL_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>. 
 */
void VECTOR_DBL_Resize(    VECTOR_DBL*    vec, 
                           int            size );

/*
 *  FUNCTION:  VECTOR_DBL_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>. 
 */
void VECTOR_DBL_GrowTo(   VECTOR_DBL*     vec, 
                           int            size );

/*
 *  FUNCTION:  VECTOR_DBL_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
void VECTOR_DBL_Push(   VECTOR_DBL*   vec, 
                        DBL           val );

/*
 *  FUNCTION:  VECTOR_DBL_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
void VECTOR_DBL_Pushback(  VECTOR_DBL*   vec, 
                           DBL           val );

/*
 *  FUNCTION:  VECTOR_DBL_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data. 
 */
DBL VECTOR_DBL_Pop( VECTOR_DBL*   vec );

/*
 *  FUNCTION:  VECTOR_DBL_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
void VECTOR_DBL_Set(   VECTOR_DBL*     vec, 
                        int            idx, 
                        DBL            val );

/*
 *  FUNCTION:  VECTOR_DBL_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
DBL* VECTOR_DBL_Get(    VECTOR_DBL*   vec, 
                        int           idx );

/*
 *  FUNCTION:  VECTOR_DBL_Get_X()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return pointer to data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 *  RETURN:    Pointer to location to <vec> idx.
 */
DBL* VECTOR_DBL_Get_X(  VECTOR_DBL*   vec, 
                        int           idx );

/*
 *  FUNCTION:  VECTOR_DBL_Get_Size()
 *  SYNOPSIS:  Get utilized length of <vec>.
 */
DBL VECTOR_DBL_Get_Size(   VECTOR_DBL*   vec );

/*
 *  FUNCTION:  VECTOR_DBL_Set_Size()
 *  SYNOPSIS:  Set utilized length of <vec>
 *  RETURN:    Pointer to location to <vec> idx.
 */
void VECTOR_DBL_Set_Size(  VECTOR_DBL*   vec, 
                           int           size );

/*
 *  FUNCTION:  VECTOR_DBL_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance found.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_DBL_Search(  VECTOR_DBL*   vec, 
                        DBL           val );

/*
 *  FUNCTION:  VECTOR_DBL_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of first occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_DBL_Search_First(  VECTOR_DBL*   vec, 
                              DBL           val );

/*
 *  FUNCTION:  VECTOR_DBL_Search_Last()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the last occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of last occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_DBL_Search_Last(   VECTOR_DBL*   vec, 
                              DBL           val );

/*
 *  FUNCTION:  VECTOR_DBL_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality, 
 *             pos if <vec_A> > <vec_B>,  
 *             neg if <vec_A> < <vec_B>.
 */
int VECTOR_DBL_Compare(    VECTOR_DBL*   vec_A, 
                           VECTOR_DBL*   vec_B );

/*
 *  FUNCTION:  VECTOR_DBL_Sort()
 *  SYNOPSIS:  Sorts <vec> data array in ascending order. In-place.
 */
void VECTOR_DBL_Sort( VECTOR_DBL*    vec );

/*
 *  FUNCTION:  VECTOR_DBL_Sort_Sub()
 *  SYNOPSIS:  Subcall - Sorts <vec> data array in ascending order on range (beg,end]. 
 *             Uses quicksort until length of subarray falls below threshold, then selection sort.
 */
void VECTOR_DBL_Sort_Sub(  VECTOR_DBL*    vec,
                           int            beg,
                           int            end );

/*
 *  FUNCTION:  VECTOR_DBL_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs selection sort on <vec> data array in ascending order on range (beg,end]. 
 */
void VECTOR_DBL_Sort_Sub_Selectsort(   VECTOR_DBL*    vec,
                                       int            beg,
                                       int            end );

/*
 *  FUNCTION:  VECTOR_DBL_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs quicksort on <vec> data array in ascending order on range (beg,end]. 
 */
void VECTOR_DBL_Sort_Sub_Quicksort( VECTOR_DBL*    vec,
                                    int            beg,
                                    int            end );

/*
 *  FUNCTION:  VECTOR_DBL_Swap()
 *  SYNOPSIS:  Swaps the values of <vec> at indexes <i> and <j>.
 *             Warning: Only checks for Out-of-Bounds when in DEBUG.
 */
void VECTOR_DBL_Swap(   VECTOR_DBL*    vec,
                        int            i,
                        int            j );

/*
 *  FUNCTION:  VECTOR_DBL_Reverse()
 *  SYNOPSIS:  Reverse the ordering of array.
 */
void VECTOR_DBL_Reverse(   VECTOR_DBL*    vec );

/*
 *  FUNCTION:  VECTOR_DBL_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
void VECTOR_DBL_Dump(   VECTOR_DBL*    vec,
                        FILE*          fp );

#endif 