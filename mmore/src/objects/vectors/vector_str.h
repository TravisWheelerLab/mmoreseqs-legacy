/*******************************************************************************
 *  FILE:      vector_str.c
 *  PURPOSE:   VECTOR_STR Object Functions
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _VECTOR_STR_H
#define _VECTOR_STR_H

/*! FUNCTION:  VECTOR_STR_Create()
 *  SYNOPSIS:  Allocates new VECTOR_STR object and returns pointer.
 */
VECTOR_STR* 
VECTOR_STR_Create();

/*! FUNCTION:  VECTOR_STR_Create_by_Size()
 *  SYNOPSIS:  Allocates new VECTOR_STR object at specific size and returns pointer.
 */
VECTOR_STR* 
VECTOR_STR_Create_by_Size( size_t    size );

/*! FUNCTION:  VECTOR_STR_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_STR and returns NULL pointer.
 */
VECTOR_STR* 
VECTOR_STR_Destroy( VECTOR_STR*   vec );

/*! FUNCTION:  VECTOR_STR_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_STR object by resetting size counter (no realloc) .
 */
STATUS_FLAG 
VECTOR_STR_Reuse( VECTOR_STR*   vec );

/*! FUNCTION:  VECTOR_STR_Fill()
 *  SYNOPSIS:  Fill VECTOR_STR object with val.
 */
STATUS_FLAG 
VECTOR_STR_Fill(  VECTOR_STR*   vec, 
                  STR           val );

/*! FUNCTION:  VECTOR_STR_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object. 
 *             Creates new VECTOR_STR for <dest> if <dest> is NULL.
 */
VECTOR_STR* 
VECTOR_STR_Copy(  VECTOR_STR*   dest, 
                  VECTOR_STR*   src );

/*! FUNCTION:  VECTOR_STR_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>. 
 */
STATUS_FLAG 
VECTOR_STR_Resize(   VECTOR_STR*    vec, 
                     size_t         size );

/*! FUNCTION:  VECTOR_STR_GetArray()
 *  SYNOPSIS:  Get <data> array from <vec>.
 */
STR* 
VECTOR_STR_GetArray(   VECTOR_STR*   vec );

/*! FUNCTION:  VECTOR_STR_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>. 
 */
STATUS_FLAG 
VECTOR_STR_GrowTo(   VECTOR_STR*    vec, 
                     size_t         size );

/*! FUNCTION:  VECTOR_STR_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
STR 
VECTOR_STR_Get(   VECTOR_STR*   vec, 
                  int           idx );

/*! FUNCTION:  VECTOR_STR_GetX()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return pointer to data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 *  RETURN:    Pointer to location to <vec> idx.
 */
STR* 
VECTOR_STR_GetX(     VECTOR_STR*   vec, 
                     int           idx );

/*! FUNCTION:  VECTOR_STR_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
STATUS_FLAG 
VECTOR_STR_Set(   VECTOR_STR*       vec, 
                  int               idx, 
                  STR               val );

/*! FUNCTION:  VECTOR_STR_Insert()
 *  SYNOPSIS:  Overwrite data from <vec> at the <idx> position in array to <val>. Deletes present value.
 */
STATUS_FLAG 
VECTOR_STR_Insert(   VECTOR_STR*   vec, 
                     int           idx, 
                     STR           val );

/*! FUNCTION:  VECTOR_STR_Delete()
 *  SYNOPSIS:  Overwrite data from <vec> at the <idx> position in array to <val>. Deletes present value.
 */
STATUS_FLAG 
VECTOR_STR_Delete(   VECTOR_STR*   vec, 
                     int           idx );

/*! FUNCTION:  VECTOR_STR_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
STATUS_FLAG 
VECTOR_STR_Push(  VECTOR_STR*   vec, 
                  STR           val );

/*! FUNCTION:  VECTOR_STR_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
STATUS_FLAG 
VECTOR_STR_Pushback(    VECTOR_STR*   vec, 
                        STR           val );

/*! FUNCTION:  VECTOR_STR_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data. 
 */
STR 
VECTOR_STR_Pop( VECTOR_STR*   vec );

/*! FUNCTION:  VECTOR_STR_Popback()
 *  SYNOPSIS:  Pop data from the end of <vec> data array and return data.
 *             Resize if array is less than half full.
 */
STR 
VECTOR_STR_Popback( VECTOR_STR*   vec );

/*! FUNCTION:  VECTOR_STR_Append()
 *  SYNOPSIS:  Push <append> data array of length <L> onto the end of <vec> data array. 
 */
STATUS_FLAG 
VECTOR_STR_Append(   VECTOR_STR*   vec, 
                     STR*          append,
                     size_t        L );

/*! FUNCTION:  VECTOR_STR_GetSize()
 *  SYNOPSIS:  Get utilized length of <vec>.
 */
int 
VECTOR_STR_GetSize(   VECTOR_STR*   vec );

/*! FUNCTION:  VECTOR_STR_SetSize()
 *  SYNOPSIS:  Set utilized length of <vec>
 *  RETURN:    Pointer to location to <vec> idx.
 */
STATUS_FLAG 
VECTOR_STR_SetSize( VECTOR_STR*   vec, 
                     size_t        size );

/*! FUNCTION:  VECTOR_STR_GetSizeAlloc()
 *  SYNOPSIS:  Get allocated length of <vec>.
 */
int 
VECTOR_STR_GetSizeAlloc(   VECTOR_STR*   vec );

/*! FUNCTION:  VECTOR_STR_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance found.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int 
VECTOR_STR_Search(   VECTOR_STR*   vec, 
                     STR           val );

/*! FUNCTION:  VECTOR_STR_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of first occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int 
VECTOR_STR_Search_First(   VECTOR_STR*   vec, 
                           STR           val );

/*! FUNCTION:  VECTOR_STR_Search_Last()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the last occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of last occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int 
VECTOR_STR_Search_Last( VECTOR_STR*   vec, 
                        STR           val );

/*! FUNCTION:  VECTOR_STR_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality, 
 *             pos if <vec_A> > <vec_B>,  
 *             neg if <vec_A> < <vec_B>.
 */
int 
VECTOR_STR_Compare(  VECTOR_STR*   vec_A, 
                     VECTOR_STR*   vec_B );

/*! FUNCTION:  VECTOR_STR_Sort()
 *  SYNOPSIS:  Sorts <vec> data array in ascending order. In-place.
 */
STATUS_FLAG 
VECTOR_STR_Sort( VECTOR_STR*    vec );

/*! FUNCTION:  VECTOR_STR_Sort_Sub()
 *  SYNOPSIS:  Subcall - Sorts <vec> data array in ascending order on range (beg,end]. 
 *             Uses quicksort until length of subarray falls below threshold, then selection sort.
 */
STATUS_FLAG 
VECTOR_STR_Sort_Sub(    VECTOR_STR*    vec,
                        int            beg,
                        int            end );

/*! FUNCTION:  VECTOR_STR_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs selection sort on <vec> data array in ascending order on range (beg,end]. 
 */
STATUS_FLAG 
VECTOR_STR_Sort_Sub_Selectsort(  VECTOR_STR*    vec,
                                 int            beg,
                                 int            end );

/*! FUNCTION:  VECTOR_STR_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs quicksort on <vec> data array in ascending order on range (beg,end]. 
 */
STATUS_FLAG 
VECTOR_STR_Sort_Sub_Quicksort(   VECTOR_STR*    vec,
                                 int            beg,
                                 int            end );

/*! FUNCTION:  VECTOR_STR_Op()
 *  SYNOPSIS:  Perform element-wise unary operation <op>(STR data) to each cell in <vec_in> and puts it in <vec_out>.
 *             Returns a pointer to <vec_out>.
 *             <vec_out> will be resized to size of <vec_in>.
 *             <vec_in> and <vec_out> can be the same vector. If <vec_out> is NULL, new vector will be created.
 */
VECTOR_STR* 
VECTOR_STR_Op(    VECTOR_STR*    vec_out,             /* input vector */
                  VECTOR_STR*    vec_in,              /* output vector (can be input vector) */
                  STR            (*op)(STR data) );   /* unary operation function */

/*! FUNCTION:  VECTOR_STR_Op()
 *  SYNOPSIS:  Perform element-wise binary operation <op>(STR data_1, STR data_2) to each cell in <vec_in_1, vec_in_2> and puts it in <vec_out>.
 *             <vec_in_1> and <vec_in_2> must be the same size.
 *             Returns a pointer to <vec_out>.
 *             <vec_out> will be resized to size of <vec_in>.
 *             <vec_in> and <vec_out> can be the same vector. If <vec_out> is NULL, new vector will be created.
 */
VECTOR_STR* 
VECTOR_STR_BinOp(    VECTOR_STR*    vec_out,                         /* output vector */ 
                     VECTOR_STR*    vec_in_1,                        /* first input vector */
                     VECTOR_STR*    vec_in_2,                        /* second input vector */
                     STR            (*op)(STR data_1, STR data_2) ); /* binary operation function */

/*! FUNCTION:  VECTOR_STR_Swap()
 *  SYNOPSIS:  Swaps the values of <vec> at indexes <i> and <j>.
 */
STATUS_FLAG 
VECTOR_STR_Swap(  VECTOR_STR*    vec,
                  int            i,
                  int            j );

/*! FUNCTION:  VECTOR_STR_Reverse()
 *  SYNOPSIS:  Reverse the ordering of array.
 */
STATUS_FLAG 
VECTOR_STR_Reverse(   VECTOR_STR*    vec );

/*! FUNCTION:  VECTOR_STR_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
STATUS_FLAG 
VECTOR_STR_Dump(  VECTOR_STR*    vec,
                  FILE*          fp );

/*! FUNCTION:  VECTOR_STR_Dump_byOpt()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer. Non-optimized.
 */
STATUS_FLAG 
VECTOR_STR_Dump_byOpt(     VECTOR_STR*    vec,
                           STR            delim,
                           STR            header,
                           FILE*          fp );

/*! FUNCTION:  VECTOR_STR_Unit_Test()
 *  SYNOPSIS:  Perform unit test for VECTOR_STR.
 */
STATUS_FLAG 
VECTOR_STR_Unit_Test();

#endif 