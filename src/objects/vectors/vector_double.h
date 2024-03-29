/*******************************************************************************
 *  - FILE:      vector_double.c
 *  - DESC:    VECTOR_DBL Object Functions
 *******************************************************************************/

#ifndef _VECTOR_DBL_H
#define _VECTOR_DBL_H

/*! FUNCTION:  VECTOR_DBL_Create()
 *  SYNOPSIS:  Allocates new VECTOR_DBL object and returns pointer.
 */
VECTOR_DBL* VECTOR_DBL_Create();

/*! FUNCTION:  VECTOR_DBL_Create_by_Size()
 *  SYNOPSIS:  Allocates new VECTOR_DBL object at specific size and returns
 * pointer.
 */
VECTOR_DBL* VECTOR_DBL_Create_by_Size(size_t size);

/*! FUNCTION:  VECTOR_DBL_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_DBL and returns NULL
 * pointer.
 */
VECTOR_DBL* VECTOR_DBL_Destroy(VECTOR_DBL* vec);

/*! FUNCTION:  VECTOR_DBL_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_DBL object by resetting size counter (no realloc) .
 */
STATUS_FLAG
VECTOR_DBL_Reuse(VECTOR_DBL* vec);

/*! FUNCTION:  VECTOR_DBL_Fill()
 *  SYNOPSIS:  Fill VECTOR_DBL object with val.
 */
STATUS_FLAG
VECTOR_DBL_Fill(VECTOR_DBL* vec, DBL val);

/*! FUNCTION:  VECTOR_DBL_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object.
 *             Creates new VECTOR_DBL for <dest> if <dest> is NULL.
 */
VECTOR_DBL* VECTOR_DBL_Copy(VECTOR_DBL* dest, VECTOR_DBL* src);

/*! FUNCTION:  VECTOR_DBL_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>.
 */
STATUS_FLAG
VECTOR_DBL_Resize(VECTOR_DBL* vec, size_t size);

/*! FUNCTION:  VECTOR_DBL_GetArray()
 *  SYNOPSIS:  Get <data> array from <vec>.
 */
DBL* VECTOR_DBL_GetArray(VECTOR_DBL* vec);

/*! FUNCTION:  VECTOR_DBL_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>.
 */
STATUS_FLAG
VECTOR_DBL_GrowTo(VECTOR_DBL* vec, size_t size);

/*! FUNCTION:  VECTOR_DBL_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return
 * data. Warning: Out-of-Bounds only checked in DEBUG.
 */
DBL VECTOR_DBL_Get(VECTOR_DBL* vec, int idx);

/*! FUNCTION:  VECTOR_DBL_GetX()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return
 * pointer to data. Warning: Out-of-Bounds only checked in DEBUG. RETURN:
 * Pointer to location to <vec> idx.
 */
DBL* VECTOR_DBL_GetX(VECTOR_DBL* vec, int idx);

/*! FUNCTION:  VECTOR_DBL_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
STATUS_FLAG
VECTOR_DBL_Set(VECTOR_DBL* vec, int idx, DBL val);

/*! FUNCTION:  VECTOR_DBL_Insert()
 *  SYNOPSIS:  Overwrite data from <vec> at the <idx> position in array to
 * <val>. Deletes present value.
 */
STATUS_FLAG
VECTOR_DBL_Insert(VECTOR_DBL* vec, int idx, DBL val);

/*! FUNCTION:  VECTOR_DBL_Delete()
 *  SYNOPSIS:  Overwrite data from <vec> at the <idx> position in array to
 * <val>. Deletes present value.
 */
STATUS_FLAG
VECTOR_DBL_Delete(VECTOR_DBL* vec, int idx);

/*! FUNCTION:  VECTOR_DBL_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
STATUS_FLAG
VECTOR_DBL_Push(VECTOR_DBL* vec, DBL val);

/*! FUNCTION:  VECTOR_DBL_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
STATUS_FLAG
VECTOR_DBL_Pushback(VECTOR_DBL* vec, DBL val);

/*! FUNCTION:  VECTOR_DBL_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data.
 */
DBL VECTOR_DBL_Pop(VECTOR_DBL* vec);

/*! FUNCTION:  VECTOR_DBL_Popback()
 *  SYNOPSIS:  Pop data from the end of <vec> data array and return data.
 *             Resize if array is less than half full.
 */
DBL VECTOR_DBL_Popback(VECTOR_DBL* vec);

/*! FUNCTION:  VECTOR_DBL_Append()
 *  SYNOPSIS:  Push <append> data array of length <L> onto the end of <vec> data
 * array.
 */
STATUS_FLAG
VECTOR_DBL_Append(VECTOR_DBL* vec, DBL* append, size_t L);

/*! FUNCTION:  VECTOR_DBL_GetSize()
 *  SYNOPSIS:  Get utilized length of <vec>.
 */
int VECTOR_DBL_GetSize(VECTOR_DBL* vec);

/*! FUNCTION:  VECTOR_DBL_SetSize()
 *  SYNOPSIS:  Set utilized length of <vec>
 *  RETURN:    Pointer to location to <vec> idx.
 */
STATUS_FLAG
VECTOR_DBL_SetSize(VECTOR_DBL* vec, size_t size);

/*! FUNCTION:  VECTOR_DBL_GetSizeAlloc()
 *  SYNOPSIS:  Get allocated length of <vec>.
 */
int VECTOR_DBL_GetSizeAlloc(VECTOR_DBL* vec);

/*! FUNCTION:  VECTOR_DBL_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array.
 *             Returns index of the first occurance found.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of occurance of <val>.
 *             Return -1 if <val> is not found.
 */
int VECTOR_DBL_Search(VECTOR_DBL* vec, DBL val);

/*! FUNCTION:  VECTOR_DBL_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array.
 *             Returns index of the first occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of first occurance of <val>.
 *             Return -1 if <val> is not found.
 */
int VECTOR_DBL_Search_First(VECTOR_DBL* vec, DBL val);

/*! FUNCTION:  VECTOR_DBL_Search_Last()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array.
 *             Returns index of the last occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of last occurance of <val>.
 *             Return -1 if <val> is not found.
 */
int VECTOR_DBL_Search_Last(VECTOR_DBL* vec, DBL val);

/*! FUNCTION:  VECTOR_DBL_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality,
 *             pos if <vec_A> > <vec_B>,
 *             neg if <vec_A> < <vec_B>.
 */
int VECTOR_DBL_Compare(VECTOR_DBL* vec_A, VECTOR_DBL* vec_B);

/*! FUNCTION:  VECTOR_DBL_Sort()
 *  SYNOPSIS:  Sorts <vec> data array in ascending order. In-place.
 */
STATUS_FLAG
VECTOR_DBL_Sort(VECTOR_DBL* vec);

/*! FUNCTION:  VECTOR_DBL_Sort_Sub()
 *  SYNOPSIS:  Subcall - Sorts <vec> data array in ascending order on range
 * (beg,end]. Uses quicksort until length of subarray falls below threshold,
 * then selection sort.
 */
STATUS_FLAG
VECTOR_DBL_Sort_Sub(VECTOR_DBL* vec, int beg, int end);

/*! FUNCTION:  VECTOR_DBL_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs selection sort on <vec> data array in ascending
 * order on range (beg,end].
 */
STATUS_FLAG
VECTOR_DBL_Sort_Sub_Selectsort(VECTOR_DBL* vec, int beg, int end);

/*! FUNCTION:  VECTOR_DBL_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs quicksort on <vec> data array in ascending order
 * on range (beg,end].
 */
STATUS_FLAG
VECTOR_DBL_Sort_Sub_Quicksort(VECTOR_DBL* vec, int beg, int end);

/*! FUNCTION:  VECTOR_DBL_Op()
 *  SYNOPSIS:  Perform element-wise unary operation <op>(DBL data) to each cell
 * in <vec_in> and puts it in <vec_out>. Returns a pointer to <vec_out>.
 *             <vec_out> will be resized to size of <vec_in>.
 *             <vec_in> and <vec_out> can be the same vector. If <vec_out> is
 * NULL, new vector will be created.
 */
VECTOR_DBL* VECTOR_DBL_Op(
    VECTOR_DBL* vec_out,  /* input vector */
    VECTOR_DBL* vec_in,   /* output vector (can be input vector) */
    DBL (*op)(DBL data)); /* unary operation function */

/*! FUNCTION:  VECTOR_DBL_Op()
 *  SYNOPSIS:  Perform element-wise binary operation <op>(DBL data_1, DBL
 * data_2) to each cell in <vec_in_1, vec_in_2> and puts it in <vec_out>.
 *             <vec_in_1> and <vec_in_2> must be the same size.
 *             Returns a pointer to <vec_out>.
 *             <vec_out> will be resized to size of <vec_in>.
 *             <vec_in> and <vec_out> can be the same vector. If <vec_out> is
 * NULL, new vector will be created.
 */
VECTOR_DBL* VECTOR_DBL_BinOp(
    VECTOR_DBL* vec_out,                /* output vector */
    VECTOR_DBL* vec_in_1,               /* first input vector */
    VECTOR_DBL* vec_in_2,               /* second input vector */
    DBL (*op)(DBL data_1, DBL data_2)); /* binary operation function */

/*! FUNCTION:  VECTOR_DBL_Swap()
 *  SYNOPSIS:  Swaps the values of <vec> at indexes <i> and <j>.
 */
STATUS_FLAG
VECTOR_DBL_Swap(VECTOR_DBL* vec, int i, int j);

/*! FUNCTION:  VECTOR_DBL_Reverse()
 *  SYNOPSIS:  Reverse the ordering of array.
 */
STATUS_FLAG
VECTOR_DBL_Reverse(VECTOR_DBL* vec);

/*! FUNCTION:  VECTOR_DBL_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
STATUS_FLAG
VECTOR_DBL_Dump(VECTOR_DBL* vec, FILE* fp);

/*! FUNCTION:  VECTOR_DBL_Dump_byOpt()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer. Non-optimized.
 */
STATUS_FLAG
VECTOR_DBL_Dump_byOpt(VECTOR_DBL* vec, STR delim, STR header, FILE* fp);

/*! FUNCTION:  VECTOR_DBL_Unit_Test()
 *  SYNOPSIS:  Perform unit test for VECTOR_DBL.
 */
STATUS_FLAG
VECTOR_DBL_Unit_Test();

#endif
