/*******************************************************************************
 *  FILE:      vector_char.c
 *  PURPOSE:   VECTOR_CHAR Object Functions
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _VECTOR_CHAR_H
#define _VECTOR_CHAR_H

/*! FUNCTION:  VECTOR_CHAR_Create()
 *  SYNOPSIS:  Allocates new VECTOR_CHAR object and returns pointer.
 */
VECTOR_CHAR* 
VECTOR_CHAR_Create();

/*! FUNCTION:  VECTOR_CHAR_Create_by_Size()
 *  SYNOPSIS:  Allocates new VECTOR_CHAR object at specific size and returns pointer.
 */
VECTOR_CHAR* 
VECTOR_CHAR_Create_by_Size( size_t    size );

/*! FUNCTION:  VECTOR_CHAR_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_CHAR and returns NULL pointer.
 */
VECTOR_CHAR* 
VECTOR_CHAR_Destroy( VECTOR_CHAR*   vec );

/*! FUNCTION:  VECTOR_CHAR_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_CHAR object by resetting size counter (no realloc) .
 */
STATUS_FLAG 
VECTOR_CHAR_Reuse( VECTOR_CHAR*   vec );

/*! FUNCTION:  VECTOR_CHAR_Fill()
 *  SYNOPSIS:  Fill VECTOR_CHAR object with val.
 */
STATUS_FLAG 
VECTOR_CHAR_Fill(  VECTOR_CHAR*   vec, 
                  CHAR           val );

/*! FUNCTION:  VECTOR_CHAR_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object. 
 *             Creates new VECTOR_CHAR for <dest> if <dest> is NULL.
 */
VECTOR_CHAR* 
VECTOR_CHAR_Copy(  VECTOR_CHAR*   dest, 
                  VECTOR_CHAR*   src );

/*! FUNCTION:  VECTOR_CHAR_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>. 
 */
STATUS_FLAG 
VECTOR_CHAR_Resize(   VECTOR_CHAR*    vec, 
                     size_t         size );

/*! FUNCTION:  VECTOR_CHAR_GetArray()
 *  SYNOPSIS:  Get <data> array from <vec>.
 */
CHAR* 
VECTOR_CHAR_GetArray(   VECTOR_CHAR*   vec );

/*! FUNCTION:  VECTOR_CHAR_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>. 
 */
STATUS_FLAG 
VECTOR_CHAR_GrowTo(   VECTOR_CHAR*    vec, 
                     size_t         size );

/*! FUNCTION:  VECTOR_CHAR_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
CHAR 
VECTOR_CHAR_Get(   VECTOR_CHAR*   vec, 
                  int           idx );

/*! FUNCTION:  VECTOR_CHAR_GetX()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return pointer to data.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 *  RETURN:    Pointer to location to <vec> idx.
 */
CHAR* 
VECTOR_CHAR_GetX(     VECTOR_CHAR*   vec, 
                     int           idx );

/*! FUNCTION:  VECTOR_CHAR_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
STATUS_FLAG 
VECTOR_CHAR_Set(   VECTOR_CHAR*       vec, 
                  int               idx, 
                  CHAR               val );

/*! FUNCTION:  VECTOR_CHAR_Insert()
 *  SYNOPSIS:  Overwrite data from <vec> at the <idx> position in array to <val>. Deletes present value.
 */
STATUS_FLAG 
VECTOR_CHAR_Insert(   VECTOR_CHAR*   vec, 
                     int           idx, 
                     CHAR           val );

/*! FUNCTION:  VECTOR_CHAR_Delete()
 *  SYNOPSIS:  Overwrite data from <vec> at the <idx> position in array to <val>. Deletes present value.
 */
STATUS_FLAG 
VECTOR_CHAR_Delete(   VECTOR_CHAR*   vec, 
                     int           idx );

/*! FUNCTION:  VECTOR_CHAR_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
STATUS_FLAG 
VECTOR_CHAR_Push(  VECTOR_CHAR*   vec, 
                  CHAR           val );

/*! FUNCTION:  VECTOR_CHAR_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
STATUS_FLAG 
VECTOR_CHAR_Pushback(    VECTOR_CHAR*   vec, 
                        CHAR           val );

/*! FUNCTION:  VECTOR_CHAR_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data. 
 */
CHAR 
VECTOR_CHAR_Pop( VECTOR_CHAR*   vec );

/*! FUNCTION:  VECTOR_CHAR_Popback()
 *  SYNOPSIS:  Pop data from the end of <vec> data array and return data.
 *             Resize if array is less than half full.
 */
CHAR 
VECTOR_CHAR_Popback( VECTOR_CHAR*   vec );

/*! FUNCTION:  VECTOR_CHAR_Append()
 *  SYNOPSIS:  Push <append> data array of length <L> onto the end of <vec> data array. 
 */
STATUS_FLAG 
VECTOR_CHAR_Append(   VECTOR_CHAR*   vec, 
                     CHAR*          append,
                     size_t        L );

/*! FUNCTION:  VECTOR_CHAR_GetSize()
 *  SYNOPSIS:  Get utilized length of <vec>.
 */
int 
VECTOR_CHAR_GetSize(   VECTOR_CHAR*   vec );

/*! FUNCTION:  VECTOR_CHAR_SetSize()
 *  SYNOPSIS:  Set utilized length of <vec>
 *  RETURN:    Pointer to location to <vec> idx.
 */
STATUS_FLAG 
VECTOR_CHAR_SetSize( VECTOR_CHAR*   vec, 
                     size_t        size );

/*! FUNCTION:  VECTOR_CHAR_GetSizeAlloc()
 *  SYNOPSIS:  Get allocated length of <vec>.
 */
int 
VECTOR_CHAR_GetSizeAlloc(   VECTOR_CHAR*   vec );

/*! FUNCTION:  VECTOR_CHAR_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance found.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int 
VECTOR_CHAR_Search(   VECTOR_CHAR*   vec, 
                     CHAR           val );

/*! FUNCTION:  VECTOR_CHAR_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the first occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of first occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int 
VECTOR_CHAR_Search_First(   VECTOR_CHAR*   vec, 
                           CHAR           val );

/*! FUNCTION:  VECTOR_CHAR_Search_Last()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array. 
 *             Returns index of the last occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of last occurance of <val>.  
 *             Return -1 if <val> is not found.
 */
int 
VECTOR_CHAR_Search_Last( VECTOR_CHAR*   vec, 
                        CHAR           val );

/*! FUNCTION:  VECTOR_CHAR_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality, 
 *             pos if <vec_A> > <vec_B>,  
 *             neg if <vec_A> < <vec_B>.
 */
int 
VECTOR_CHAR_Compare(  VECTOR_CHAR*   vec_A, 
                     VECTOR_CHAR*   vec_B );

/*! FUNCTION:  VECTOR_CHAR_Sort()
 *  SYNOPSIS:  Sorts <vec> data array in ascending order. In-place.
 */
STATUS_FLAG 
VECTOR_CHAR_Sort( VECTOR_CHAR*    vec );

/*! FUNCTION:  VECTOR_CHAR_Sort_Sub()
 *  SYNOPSIS:  Subcall - Sorts <vec> data array in ascending order on range (beg,end]. 
 *             Uses quicksort until length of subarray falls below threshold, then selection sort.
 */
STATUS_FLAG 
VECTOR_CHAR_Sort_Sub(    VECTOR_CHAR*    vec,
                        int            beg,
                        int            end );

/*! FUNCTION:  VECTOR_CHAR_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs selection sort on <vec> data array in ascending order on range (beg,end]. 
 */
STATUS_FLAG 
VECTOR_CHAR_Sort_Sub_Selectsort(  VECTOR_CHAR*    vec,
                                 int            beg,
                                 int            end );

/*! FUNCTION:  VECTOR_CHAR_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs quicksort on <vec> data array in ascending order on range (beg,end]. 
 */
STATUS_FLAG 
VECTOR_CHAR_Sort_Sub_Quicksort(   VECTOR_CHAR*    vec,
                                 int            beg,
                                 int            end );

/*! FUNCTION:  VECTOR_CHAR_Op()
 *  SYNOPSIS:  Perform element-wise unary operation <op>(CHAR data) to each cell in <vec_in> and puts it in <vec_out>.
 *             Returns a pointer to <vec_out>.
 *             <vec_out> will be resized to size of <vec_in>.
 *             <vec_in> and <vec_out> can be the same vector. If <vec_out> is NULL, new vector will be created.
 */
VECTOR_CHAR* 
VECTOR_CHAR_Op(    VECTOR_CHAR*    vec_out,             /* input vector */
                  VECTOR_CHAR*    vec_in,              /* output vector (can be input vector) */
                  CHAR            (*op)(CHAR data) );   /* unary operation function */

/*! FUNCTION:  VECTOR_CHAR_Op()
 *  SYNOPSIS:  Perform element-wise binary operation <op>(CHAR data_1, CHAR data_2) to each cell in <vec_in_1, vec_in_2> and puts it in <vec_out>.
 *             <vec_in_1> and <vec_in_2> must be the same size.
 *             Returns a pointer to <vec_out>.
 *             <vec_out> will be resized to size of <vec_in>.
 *             <vec_in> and <vec_out> can be the same vector. If <vec_out> is NULL, new vector will be created.
 */
VECTOR_CHAR* 
VECTOR_CHAR_BinOp(    VECTOR_CHAR*    vec_out,                         /* output vector */ 
                     VECTOR_CHAR*    vec_in_1,                        /* first input vector */
                     VECTOR_CHAR*    vec_in_2,                        /* second input vector */
                     CHAR            (*op)(CHAR data_1, CHAR data_2) ); /* binary operation function */

/*! FUNCTION:  VECTOR_CHAR_Swap()
 *  SYNOPSIS:  Swaps the values of <vec> at indexes <i> and <j>.
 */
STATUS_FLAG 
VECTOR_CHAR_Swap(  VECTOR_CHAR*    vec,
                  int            i,
                  int            j );

/*! FUNCTION:  VECTOR_CHAR_Reverse()
 *  SYNOPSIS:  Reverse the ordering of array.
 */
STATUS_FLAG 
VECTOR_CHAR_Reverse(   VECTOR_CHAR*    vec );

/*! FUNCTION:  VECTOR_CHAR_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
STATUS_FLAG 
VECTOR_CHAR_Dump(  VECTOR_CHAR*    vec,
                  FILE*          fp );

/*! FUNCTION:  VECTOR_CHAR_Dump_byOpt()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer. Non-optimized.
 */
STATUS_FLAG 
VECTOR_CHAR_Dump_byOpt(     VECTOR_CHAR*    vec,
                           STR            delim,
                           STR            header,
                           FILE*          fp );

/*! FUNCTION:  VECTOR_CHAR_Unit_Test()
 *  SYNOPSIS:  Perform unit test for VECTOR_CHAR.
 */
STATUS_FLAG 
VECTOR_CHAR_Unit_Test();

#endif 