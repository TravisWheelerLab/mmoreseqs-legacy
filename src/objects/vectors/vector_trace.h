/*******************************************************************************
 *  - FILE:   vector_trace.c
 *  - DESC:    VECTOR_TRACE Object Functions
 *******************************************************************************/

#ifndef _VECTOR_TRACE_H
#define _VECTOR_TRACE_H

/*! FUNCTION:  VECTOR_TRACE_Create()
 *  SYNOPSIS:  Allocates new VECTOR_TRACE object and returns pointer.
 */
VECTOR_TRACE* VECTOR_TRACE_Create();

/*! FUNCTION:  VECTOR_TRACE_Create_by_Size()
 *  SYNOPSIS:  Allocates new VECTOR_TRACE object at specific size and returns
 * pointer.
 */
VECTOR_TRACE* VECTOR_TRACE_Create_by_Size(size_t size);

/*! FUNCTION:  VECTOR_TRACE_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_TRACE and returns NULL
 * pointer.
 */
VECTOR_TRACE* VECTOR_TRACE_Destroy(VECTOR_TRACE* vec);

/*! FUNCTION:  VECTOR_TRACE_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_TRACE object by resetting size counter (no realloc)
 * .
 */
STATUS_FLAG
VECTOR_TRACE_Reuse(VECTOR_TRACE* vec);

/*! FUNCTION:  VECTOR_TRACE_Fill()
 *  SYNOPSIS:  Fill VECTOR_TRACE object with val.
 */
STATUS_FLAG
VECTOR_TRACE_Fill(VECTOR_TRACE* vec, TRACE val);

/*! FUNCTION:  VECTOR_TRACE_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object.
 *             Creates new VECTOR_TRACE for <dest> if <dest> is NULL.
 */
VECTOR_TRACE* VECTOR_TRACE_Copy(VECTOR_TRACE* dest, VECTOR_TRACE* src);

/*! FUNCTION:  VECTOR_TRACE_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>.
 */
STATUS_FLAG
VECTOR_TRACE_Resize(VECTOR_TRACE* vec, size_t size);

/*! FUNCTION:  VECTOR_TRACE_GetArray()
 *  SYNOPSIS:  Get <data> array from <vec>.
 */
TRACE* VECTOR_TRACE_GetArray(VECTOR_TRACE* vec);

/*! FUNCTION:  VECTOR_TRACE_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>.
 */
STATUS_FLAG
VECTOR_TRACE_GrowTo(VECTOR_TRACE* vec, size_t size);

/*! FUNCTION:  VECTOR_TRACE_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return
 * data. Warning: Out-of-Bounds only checked in DEBUG.
 */
TRACE
VECTOR_TRACE_Get(VECTOR_TRACE* vec, int idx);

/*! FUNCTION:  VECTOR_TRACE_GetX()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return
 * pointer to data. Warning: Out-of-Bounds only checked in DEBUG. RETURN:
 * Pointer to location to <vec> idx.
 */
TRACE* VECTOR_TRACE_GetX(VECTOR_TRACE* vec, int idx);

/*! FUNCTION:  VECTOR_TRACE_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
STATUS_FLAG
VECTOR_TRACE_Set(VECTOR_TRACE* vec, int idx, TRACE val);

/*! FUNCTION:  VECTOR_TRACE_Insert()
 *  SYNOPSIS:  Overwrite data from <vec> at the <idx> position in array to
 * <val>. Deletes present value.
 */
STATUS_FLAG
VECTOR_TRACE_Insert(VECTOR_TRACE* vec, int idx, TRACE val);

/*! FUNCTION:  VECTOR_TRACE_Delete()
 *  SYNOPSIS:  Overwrite data from <vec> at the <idx> position in array to
 * <val>. Deletes present value.
 */
STATUS_FLAG
VECTOR_TRACE_Delete(VECTOR_TRACE* vec, int idx);

/*! FUNCTION:  VECTOR_TRACE_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
STATUS_FLAG
VECTOR_TRACE_Push(VECTOR_TRACE* vec, TRACE val);

/*! FUNCTION:  VECTOR_TRACE_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
STATUS_FLAG
VECTOR_TRACE_Pushback(VECTOR_TRACE* vec, TRACE val);

/*! FUNCTION:  VECTOR_TRACE_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data.
 */
TRACE
VECTOR_TRACE_Pop(VECTOR_TRACE* vec);

/*! FUNCTION:  VECTOR_TRACE_Popback()
 *  SYNOPSIS:  Pop data from the end of <vec> data array and return data.
 *             Resize if array is less than half full.
 */
TRACE
VECTOR_TRACE_Popback(VECTOR_TRACE* vec);

/*! FUNCTION:  VECTOR_TRACE_Append()
 *  SYNOPSIS:  Push <append> data array of length <L> onto the end of <vec> data
 * array.
 */
STATUS_FLAG
VECTOR_TRACE_Append(VECTOR_TRACE* vec, TRACE* append, size_t L);

/*! FUNCTION:  VECTOR_TRACE_GetSize()
 *  SYNOPSIS:  Get utilized length of <vec>.
 */
int VECTOR_TRACE_GetSize(VECTOR_TRACE* vec);

/*! FUNCTION:  VECTOR_TRACE_SetSize()
 *  SYNOPSIS:  Set utilized length of <vec>
 *  RETURN:    Pointer to location to <vec> idx.
 */
STATUS_FLAG
VECTOR_TRACE_SetSize(VECTOR_TRACE* vec, size_t size);

/*! FUNCTION:  VECTOR_TRACE_GetSizeAlloc()
 *  SYNOPSIS:  Get allocated length of <vec>.
 */
int VECTOR_TRACE_GetSizeAlloc(VECTOR_TRACE* vec);

/*! FUNCTION:  VECTOR_TRACE_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array.
 *             Returns index of the first occurance found.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of occurance of <val>.
 *             Return -1 if <val> is not found.
 */
int VECTOR_TRACE_Search(VECTOR_TRACE* vec, TRACE val);

/*! FUNCTION:  VECTOR_TRACE_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array.
 *             Returns index of the first occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of first occurance of <val>.
 *             Return -1 if <val> is not found.
 */
int VECTOR_TRACE_Search_First(VECTOR_TRACE* vec, TRACE val);

/*! FUNCTION:  VECTOR_TRACE_Search_Last()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array.
 *             Returns index of the last occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of last occurance of <val>.
 *             Return -1 if <val> is not found.
 */
int VECTOR_TRACE_Search_Last(VECTOR_TRACE* vec, TRACE val);

/*! FUNCTION:  VECTOR_TRACE_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality,
 *             pos if <vec_A> > <vec_B>,
 *             neg if <vec_A> < <vec_B>.
 */
int VECTOR_TRACE_Compare(VECTOR_TRACE* vec_A, VECTOR_TRACE* vec_B);

/*! FUNCTION:  VECTOR_TRACE_Sort()
 *  SYNOPSIS:  Sorts <vec> data array in ascending order. In-place.
 */
STATUS_FLAG
VECTOR_TRACE_Sort(VECTOR_TRACE* vec);

/*! FUNCTION:  VECTOR_TRACE_Sort_Sub()
 *  SYNOPSIS:  Subcall - Sorts <vec> data array in ascending order on range
 * (beg,end]. Uses quicksort until length of subarray falls below threshold,
 * then selection sort.
 */
STATUS_FLAG
VECTOR_TRACE_Sort_Sub(VECTOR_TRACE* vec, int beg, int end);

/*! FUNCTION:  VECTOR_TRACE_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs selection sort on <vec> data array in ascending
 * order on range (beg,end].
 */
STATUS_FLAG
VECTOR_TRACE_Sort_Sub_Selectsort(VECTOR_TRACE* vec, int beg, int end);

/*! FUNCTION:  VECTOR_TRACE_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs quicksort on <vec> data array in ascending order
 * on range (beg,end].
 */
STATUS_FLAG
VECTOR_TRACE_Sort_Sub_Quicksort(VECTOR_TRACE* vec, int beg, int end);

/*! FUNCTION:  VECTOR_TRACE_Op()
 *  SYNOPSIS:  Perform element-wise unary operation <op>(TRACE data) to each
 * cell in <vec_in> and puts it in <vec_out>. Returns a pointer to <vec_out>.
 *             <vec_out> will be resized to size of <vec_in>.
 *             <vec_in> and <vec_out> can be the same vector. If <vec_out> is
 * NULL, new vector will be created.
 */
VECTOR_TRACE* VECTOR_TRACE_Op(
    VECTOR_TRACE* vec_out,    /* input vector */
    VECTOR_TRACE* vec_in,     /* output vector (can be input vector) */
    TRACE (*op)(TRACE data)); /* unary operation function */

/*! FUNCTION:  VECTOR_TRACE_Op()
 *  SYNOPSIS:  Perform element-wise binary operation <op>(TRACE data_1, TRACE
 * data_2) to each cell in <vec_in_1, vec_in_2> and puts it in <vec_out>.
 *             <vec_in_1> and <vec_in_2> must be the same size.
 *             Returns a pointer to <vec_out>.
 *             <vec_out> will be resized to size of <vec_in>.
 *             <vec_in> and <vec_out> can be the same vector. If <vec_out> is
 * NULL, new vector will be created.
 */
VECTOR_TRACE* VECTOR_TRACE_BinOp(
    VECTOR_TRACE* vec_out,                    /* output vector */
    VECTOR_TRACE* vec_in_1,                   /* first input vector */
    VECTOR_TRACE* vec_in_2,                   /* second input vector */
    TRACE (*op)(TRACE data_1, TRACE data_2)); /* binary operation function */

/*! FUNCTION:  VECTOR_TRACE_Swap()
 *  SYNOPSIS:  Swaps the values of <vec> at indexes <i> and <j>.
 */
STATUS_FLAG
VECTOR_TRACE_Swap(VECTOR_TRACE* vec, int i, int j);

/*! FUNCTION:  VECTOR_TRACE_Reverse()
 *  SYNOPSIS:  Reverse the ordering of array.
 */
STATUS_FLAG
VECTOR_TRACE_Reverse(VECTOR_TRACE* vec);

/*! FUNCTION:  VECTOR_TRACE_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
STATUS_FLAG
VECTOR_TRACE_Dump(VECTOR_TRACE* vec, FILE* fp);

/*! FUNCTION:  VECTOR_TRACE_Dump_byOpt()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer. Non-optimized.
 */
STATUS_FLAG
VECTOR_TRACE_Dump_byOpt(VECTOR_TRACE* vec, STR delim, STR header, FILE* fp);

/*! FUNCTION:  VECTOR_TRACE_Unit_Test()
 *  SYNOPSIS:  Perform unit test for VECTOR_TRACE.
 */
STATUS_FLAG
VECTOR_TRACE_Unit_Test();

#endif
