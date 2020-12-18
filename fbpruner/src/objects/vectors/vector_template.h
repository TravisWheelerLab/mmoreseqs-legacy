/*******************************************************************************
 *  FILE:      vector_template.c
 *  PURPOSE:   VECTOR_XXX Object Functions
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _VECTOR_XXX_H
#define _VECTOR_XXX_H

/*
 *  FUNCTION:  VECTOR_XXX_Create()
 *  SYNOPSIS:  Allocates new VECTOR_XXX object and returns pointer.
 */
VECTOR_XXX* VECTOR_XXX_Create();

/*
 *  FUNCTION:  VECTOR_XXX_Create()
 *  SYNOPSIS:  Allocates new VECTOR_XXX object at specific size and returns pointer.
 */
VECTOR_XXX* VECTOR_XXX_Create_by_Size( int    size );

/*
 *  FUNCTION:  VECTOR_XXX_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_XXX.
 */
void* VECTOR_XXX_Destroy( VECTOR_XXX*   vec );

/*
 *  FUNCTION:  VECTOR_XXX_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_XXX object by resetting size counter (no realloc) .
 */
void VECTOR_XXX_Reuse( VECTOR_XXX*   vec );

/*
 *  FUNCTION:  VECTOR_XXX_Fill()
 *  SYNOPSIS:  Fill VECTOR_XXX object with val.
 */
void VECTOR_XXX_Fill(   VECTOR_XXX*   vec, 
                        XXX           val );

/*
 *  FUNCTION:  VECTOR_XXX_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object. 
 *             Creates new VECTOR_XXX for <dest> if <dest> is NULL.
 */
VECTOR_XXX* VECTOR_XXX_Copy(  VECTOR_XXX*   src, 
                              VECTOR_XXX*   dest );

/*
 *  FUNCTION:  VECTOR_XXX_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>. 
 */
void VECTOR_XXX_Resize(    VECTOR_XXX*    vec, 
                           int            size );

/*
 *  FUNCTION:  VECTOR_XXX_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>. 
 */
void VECTOR_XXX_GrowTo(   VECTOR_XXX*     vec, 
                           int            size );

/*
 *  FUNCTION:  VECTOR_XXX_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
void VECTOR_XXX_Push(   VECTOR_XXX*   vec, 
                        XXX           val );

/*
 *  FUNCTION:  VECTOR_XXX_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
void VECTOR_XXX_Pushback(  VECTOR_XXX*   vec, 
                           XXX           val );

/*
 *  FUNCTION:  VECTOR_XXX_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data. 
 */
XXX VECTOR_XXX_Pop( VECTOR_XXX*   vec );

/*
 *  FUNCTION:  VECTOR_XXX_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
void VECTOR_XXX_Set(   VECTOR_XXX*     vec, 
                        int            idx, 
                        XXX            val );

/*
 *  FUNCTION:  VECTOR_XXX_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
XXX VECTOR_XXX_Get(    VECTOR_XXX*   vec, 
                        int           idx );

/*
 *  FUNCTION:  VECTOR_XXX_Get_X()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return pointer to data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 *  RETURN:    Pointer to location to <vec> idx.
 */
XXX* VECTOR_XXX_Get_X(  VECTOR_XXX*   vec, 
                        int           idx );

/*
 *  FUNCTION:  VECTOR_XXX_Get_Size()
 *  SYNOPSIS:  Get utilized length of <vec>.
 */
int VECTOR_XXX_Get_Size(   VECTOR_XXX*   vec );

/*
 *  FUNCTION:  VECTOR_XXX_Set_Size()
 *  SYNOPSIS:  Set utilized length of <vec>
 *  RETURN:    Pointer to location to <vec> idx.
 */
void VECTOR_XXX_Set_Size(  VECTOR_XXX*   vec, 
                           int           size );

/*
 *  FUNCTION:  VECTOR_XXX_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance found.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_XXX_Search(  VECTOR_XXX*   vec, 
                        XXX           val );

/*
 *  FUNCTION:  VECTOR_XXX_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of first occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_XXX_Search_First(  VECTOR_XXX*   vec, 
                              XXX           val );

/*
 *  FUNCTION:  VECTOR_XXX_Search_Last()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the last occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of last occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int VECTOR_XXX_Search_Last(   VECTOR_XXX*   vec, 
                              XXX           val );

/*
 *  FUNCTION:  VECTOR_XXX_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality, 
 *             pos if <vec_A> > <vec_B>,  
 *             neg if <vec_A> < <vec_B>.
 */
int VECTOR_XXX_Compare(    VECTOR_XXX*   vec_A, 
                           VECTOR_XXX*   vec_B );

/*
 *  FUNCTION:  VECTOR_XXX_Sort()
 *  SYNOPSIS:  Sorts <vec> data array in ascending order. In-place.
 */
void VECTOR_XXX_Sort( VECTOR_XXX*    vec );

/*
 *  FUNCTION:  VECTOR_XXX_Sort_Sub()
 *  SYNOPSIS:  Subcall - Sorts <vec> data array in ascending order on range (beg,end]. 
 *             Uses quicksort until length of subarray falls below threshold, then selection sort.
 */
void VECTOR_XXX_Sort_Sub(  VECTOR_XXX*    vec,
                           int            beg,
                           int            end );

/*
 *  FUNCTION:  VECTOR_XXX_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs selection sort on <vec> data array in ascending order on range (beg,end]. 
 */
void VECTOR_XXX_Sort_Sub_Selectsort(   VECTOR_XXX*    vec,
                                       int            beg,
                                       int            end );

/*
 *  FUNCTION:  VECTOR_XXX_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs quicksort on <vec> data array in ascending order on range (beg,end]. 
 */
void VECTOR_XXX_Sort_Sub_Quicksort( VECTOR_XXX*    vec,
                                    int            beg,
                                    int            end );

/*
 *  FUNCTION:  VECTOR_XXX_Swap()
 *  SYNOPSIS:  Swaps the values of <vec> at indexes <i> and <j>.
 *             Warning: Only checks for Out-of-Bounds when in DEBUG.
 */
void VECTOR_XXX_Swap(   VECTOR_XXX*    vec,
                        int            i,
                        int            j );

/*
 *  FUNCTION:  VECTOR_XXX_Reverse()
 *  SYNOPSIS:  Reverse the ordering of array.
 */
void VECTOR_XXX_Reverse(   VECTOR_XXX*    vec );

/*
 *  FUNCTION:  VECTOR_XXX_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
void VECTOR_XXX_Dump(   VECTOR_XXX*    vec,
                        FILE*          fp );

/*
 *  FUNCTION:  VECTOR_XXX_Unit_Test()
 *  SYNOPSIS:  Perform unit test for VECTOR_XXX.
 */
void VECTOR_XXX_Unit_Test();

#endif 