/*******************************************************************************
 *  - FILE:   vector_int.c
 *  - DESC:    VECTOR_INT Object Functions
 *******************************************************************************/

#ifndef _VECTOR_INT_H
#define _VECTOR_INT_H

/*! FUNCTION:  VECTOR_INT_Create()
 *  SYNOPSIS:  Allocates new VECTOR_INT object and returns pointer.
 */
VECTOR_INT* VECTOR_INT_Create();

/*! FUNCTION:  VECTOR_INT_Create_by_Size()
 *  SYNOPSIS:  Allocates new VECTOR_INT object at specific size and returns
 * pointer.
 */
VECTOR_INT* VECTOR_INT_Create_by_Size(size_t size);

/*! FUNCTION:  VECTOR_INT_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_INT and returns NULL
 * pointer.
 */
VECTOR_INT* VECTOR_INT_Destroy(VECTOR_INT* vec);

/*! FUNCTION:  VECTOR_INT_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_INT object by resetting size counter (no realloc) .
 */
STATUS_FLAG
VECTOR_INT_Reuse(VECTOR_INT* vec);

/*! FUNCTION:  VECTOR_INT_Fill()
 *  SYNOPSIS:  Fill VECTOR_INT object with val.
 */
STATUS_FLAG
VECTOR_INT_Fill(VECTOR_INT* vec, INT val);

/*! FUNCTION:  VECTOR_INT_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object.
 *             Creates new VECTOR_INT for <dest> if <dest> is NULL.
 */
VECTOR_INT* VECTOR_INT_Copy(VECTOR_INT* dest, VECTOR_INT* src);

/*! FUNCTION:  VECTOR_INT_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>.
 */
STATUS_FLAG
VECTOR_INT_Resize(VECTOR_INT* vec, size_t size);

/*! FUNCTION:  VECTOR_INT_GetArray()
 *  SYNOPSIS:  Get <data> array from <vec>.
 */
INT* VECTOR_INT_GetArray(VECTOR_INT* vec);

/*! FUNCTION:  VECTOR_INT_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>.
 */
STATUS_FLAG
VECTOR_INT_GrowTo(VECTOR_INT* vec, size_t size);

/*! FUNCTION:  VECTOR_INT_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return
 * data. Warning: Out-of-Bounds only checked in DEBUG.
 */
INT VECTOR_INT_Get(VECTOR_INT* vec, int idx);

/*! FUNCTION:  VECTOR_INT_GetX()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return
 * pointer to data. Warning: Out-of-Bounds only checked in DEBUG. RETURN:
 * Pointer to location to <vec> idx.
 */
INT* VECTOR_INT_GetX(VECTOR_INT* vec, int idx);

/*! FUNCTION:  VECTOR_INT_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
STATUS_FLAG
VECTOR_INT_Set(VECTOR_INT* vec, int idx, INT val);

/*! FUNCTION:  VECTOR_INT_Insert()
 *  SYNOPSIS:  Overwrite data from <vec> at the <idx> position in array to
 * <val>. Deletes present value.
 */
STATUS_FLAG
VECTOR_INT_Insert(VECTOR_INT* vec, int idx, INT val);

/*! FUNCTION:  VECTOR_INT_Delete()
 *  SYNOPSIS:  Overwrite data from <vec> at the <idx> position in array to
 * <val>. Deletes present value.
 */
STATUS_FLAG
VECTOR_INT_Delete(VECTOR_INT* vec, int idx);

/*! FUNCTION:  VECTOR_INT_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
STATUS_FLAG
VECTOR_INT_Push(VECTOR_INT* vec, INT val);

/*! FUNCTION:  VECTOR_INT_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
STATUS_FLAG
VECTOR_INT_Pushback(VECTOR_INT* vec, INT val);

/*! FUNCTION:  VECTOR_INT_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data.
 */
INT VECTOR_INT_Pop(VECTOR_INT* vec);

/*! FUNCTION:  VECTOR_INT_Popback()
 *  SYNOPSIS:  Pop data from the end of <vec> data array and return data.
 *             Resize if array is less than half full.
 */
INT VECTOR_INT_Popback(VECTOR_INT* vec);

/*! FUNCTION:  VECTOR_INT_Append()
 *  SYNOPSIS:  Push <append> data array of length <L> onto the end of <vec> data
 * array.
 */
STATUS_FLAG
VECTOR_INT_Append(VECTOR_INT* vec, INT* append, size_t L);

/*! FUNCTION:  VECTOR_INT_GetSize()
 *  SYNOPSIS:  Get utilized length of <vec>.
 */
int VECTOR_INT_GetSize(VECTOR_INT* vec);

/*! FUNCTION:  VECTOR_INT_SetSize()
 *  SYNOPSIS:  Set utilized length of <vec>
 *  RETURN:    Pointer to location to <vec> idx.
 */
STATUS_FLAG
VECTOR_INT_SetSize(VECTOR_INT* vec, size_t size);

/*! FUNCTION:  VECTOR_INT_GetSizeAlloc()
 *  SYNOPSIS:  Get allocated length of <vec>.
 */
int VECTOR_INT_GetSizeAlloc(VECTOR_INT* vec);

/*! FUNCTION:  VECTOR_INT_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array.
 *             Returns index of the first occurance found.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of occurance of <val>.
 *             Return -1 if <val> is not found.
 */
int VECTOR_INT_Search(VECTOR_INT* vec, INT val);

/*! FUNCTION:  VECTOR_INT_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array.
 *             Returns index of the first occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of first occurance of <val>.
 *             Return -1 if <val> is not found.
 */
int VECTOR_INT_Search_First(VECTOR_INT* vec, INT val);

/*! FUNCTION:  VECTOR_INT_Search_Last()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array.
 *             Returns index of the last occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of last occurance of <val>.
 *             Return -1 if <val> is not found.
 */
int VECTOR_INT_Search_Last(VECTOR_INT* vec, INT val);

/*! FUNCTION:  VECTOR_INT_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality,
 *             pos if <vec_A> > <vec_B>,
 *             neg if <vec_A> < <vec_B>.
 */
int VECTOR_INT_Compare(VECTOR_INT* vec_A, VECTOR_INT* vec_B);

/*! FUNCTION:  VECTOR_INT_Sort()
 *  SYNOPSIS:  Sorts <vec> data array in ascending order. In-place.
 */
STATUS_FLAG
VECTOR_INT_Sort(VECTOR_INT* vec);

/*! FUNCTION:  VECTOR_INT_Sort_Sub()
 *  SYNOPSIS:  Subcall - Sorts <vec> data array in ascending order on range
 * (beg,end]. Uses quicksort until length of subarray falls below threshold,
 * then selection sort.
 */
STATUS_FLAG
VECTOR_INT_Sort_Sub(VECTOR_INT* vec, int beg, int end);

/*! FUNCTION:  VECTOR_INT_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs selection sort on <vec> data array in ascending
 * order on range (beg,end].
 */
STATUS_FLAG
VECTOR_INT_Sort_Sub_Selectsort(VECTOR_INT* vec, int beg, int end);

/*! FUNCTION:  VECTOR_INT_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs quicksort on <vec> data array in ascending order
 * on range (beg,end].
 */
STATUS_FLAG
VECTOR_INT_Sort_Sub_Quicksort(VECTOR_INT* vec, int beg, int end);

/*! FUNCTION:  VECTOR_INT_Op()
 *  SYNOPSIS:  Perform element-wise unary operation <op>(INT data) to each cell
 * in <vec_in> and puts it in <vec_out>. Returns a pointer to <vec_out>.
 *             <vec_out> will be resized to size of <vec_in>.
 *             <vec_in> and <vec_out> can be the same vector. If <vec_out> is
 * NULL, new vector will be created.
 */
VECTOR_INT* VECTOR_INT_Op(
    VECTOR_INT* vec_out,  /* input vector */
    VECTOR_INT* vec_in,   /* output vector (can be input vector) */
    INT (*op)(INT data)); /* unary operation function */

/*! FUNCTION:  VECTOR_INT_Op()
 *  SYNOPSIS:  Perform element-wise binary operation <op>(INT data_1, INT
 * data_2) to each cell in <vec_in_1, vec_in_2> and puts it in <vec_out>.
 *             <vec_in_1> and <vec_in_2> must be the same size.
 *             Returns a pointer to <vec_out>.
 *             <vec_out> will be resized to size of <vec_in>.
 *             <vec_in> and <vec_out> can be the same vector. If <vec_out> is
 * NULL, new vector will be created.
 */
VECTOR_INT* VECTOR_INT_BinOp(
    VECTOR_INT* vec_out,                /* output vector */
    VECTOR_INT* vec_in_1,               /* first input vector */
    VECTOR_INT* vec_in_2,               /* second input vector */
    INT (*op)(INT data_1, INT data_2)); /* binary operation function */

/*! FUNCTION:  VECTOR_INT_Swap()
 *  SYNOPSIS:  Swaps the values of <vec> at indexes <i> and <j>.
 */
STATUS_FLAG
VECTOR_INT_Swap(VECTOR_INT* vec, int i, int j);

/*! FUNCTION:  VECTOR_INT_Reverse()
 *  SYNOPSIS:  Reverse the ordering of array.
 */
STATUS_FLAG
VECTOR_INT_Reverse(VECTOR_INT* vec);

/*! FUNCTION:  VECTOR_INT_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
STATUS_FLAG
VECTOR_INT_Dump(VECTOR_INT* vec, FILE* fp);

/*! FUNCTION:  VECTOR_INT_Dump_byOpt()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer. Non-optimized.
 */
STATUS_FLAG
VECTOR_INT_Dump_byOpt(VECTOR_INT* vec, STR delim, STR header, FILE* fp);

/*! FUNCTION:  VECTOR_INT_Unit_Test()
 *  SYNOPSIS:  Perform unit test for VECTOR_INT.
 */
STATUS_FLAG
VECTOR_INT_Unit_Test();

#endif
