/*******************************************************************************
 *  - FILE:   vector_ptr.c
 *  - DESC:    VECTOR_PTR Object Functions
 *******************************************************************************/

#ifndef _VECTOR_PTR_H
#define _VECTOR_PTR_H

/*! FUNCTION:  VECTOR_PTR_Create()
 *  SYNOPSIS:  Allocates new VECTOR_PTR object and returns pointer.
 */
VECTOR_PTR* VECTOR_PTR_Create();

/*! FUNCTION:  VECTOR_PTR_Create_by_Size()
 *  SYNOPSIS:  Allocates new VECTOR_PTR object at specific size and returns
 * pointer.
 */
VECTOR_PTR* VECTOR_PTR_Create_by_Size(size_t size);

/*! FUNCTION:  VECTOR_PTR_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_PTR and returns NULL
 * pointer.
 */
VECTOR_PTR* VECTOR_PTR_Destroy(VECTOR_PTR* vec);

/*! FUNCTION:  VECTOR_PTR_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_PTR object by resetting size counter (no realloc) .
 */
STATUS_FLAG
VECTOR_PTR_Reuse(VECTOR_PTR* vec);

/*! FUNCTION:  VECTOR_PTR_Fill()
 *  SYNOPSIS:  Fill VECTOR_PTR object with val.
 */
STATUS_FLAG
VECTOR_PTR_Fill(VECTOR_PTR* vec, PTR val);

/*! FUNCTION:  VECTOR_PTR_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object.
 *             Creates new VECTOR_PTR for <dest> if <dest> is NULL.
 */
VECTOR_PTR* VECTOR_PTR_Copy(VECTOR_PTR* dest, VECTOR_PTR* src);

/*! FUNCTION:  VECTOR_PTR_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>.
 */
STATUS_FLAG
VECTOR_PTR_Resize(VECTOR_PTR* vec, size_t size);

/*! FUNCTION:  VECTOR_PTR_GetArray()
 *  SYNOPSIS:  Get <data> array from <vec>.
 */
PTR* VECTOR_PTR_GetArray(VECTOR_PTR* vec);

/*! FUNCTION:  VECTOR_PTR_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>.
 */
STATUS_FLAG
VECTOR_PTR_GrowTo(VECTOR_PTR* vec, size_t size);

/*! FUNCTION:  VECTOR_PTR_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return
 * data. Warning: Out-of-Bounds only checked in DEBUG.
 */
PTR VECTOR_PTR_Get(VECTOR_PTR* vec, int idx);

/*! FUNCTION:  VECTOR_PTR_GetX()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return
 * pointer to data. Warning: Out-of-Bounds only checked in DEBUG. RETURN:
 * Pointer to location to <vec> idx.
 */
PTR* VECTOR_PTR_GetX(VECTOR_PTR* vec, int idx);

/*! FUNCTION:  VECTOR_PTR_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
STATUS_FLAG
VECTOR_PTR_Set(VECTOR_PTR* vec, int idx, PTR val);

/*! FUNCTION:  VECTOR_PTR_Insert()
 *  SYNOPSIS:  Overwrite data from <vec> at the <idx> position in array to
 * <val>. Deletes present value.
 */
STATUS_FLAG
VECTOR_PTR_Insert(VECTOR_PTR* vec, int idx, PTR val);

/*! FUNCTION:  VECTOR_PTR_Delete()
 *  SYNOPSIS:  Overwrite data from <vec> at the <idx> position in array to
 * <val>. Deletes present value.
 */
STATUS_FLAG
VECTOR_PTR_Delete(VECTOR_PTR* vec, int idx);

/*! FUNCTION:  VECTOR_PTR_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
STATUS_FLAG
VECTOR_PTR_Push(VECTOR_PTR* vec, PTR val);

/*! FUNCTION:  VECTOR_PTR_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             and resize array if array is full.
 */
STATUS_FLAG
VECTOR_PTR_Pushback(VECTOR_PTR* vec, PTR val);

/*! FUNCTION:  VECTOR_PTR_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, and return data.
 */
PTR VECTOR_PTR_Pop(VECTOR_PTR* vec);

/*! FUNCTION:  VECTOR_PTR_Popback()
 *  SYNOPSIS:  Pop data from the end of <vec> data array and return data.
 *             Resize if array is less than half full.
 */
PTR VECTOR_PTR_Popback(VECTOR_PTR* vec);

/*! FUNCTION:  VECTOR_PTR_Append()
 *  SYNOPSIS:  Push <append> data array of length <L> onto the end of <vec> data
 * array.
 */
STATUS_FLAG
VECTOR_PTR_Append(VECTOR_PTR* vec, PTR* append, size_t L);

/*! FUNCTION:  VECTOR_PTR_GetSize()
 *  SYNOPSIS:  Get utilized length of <vec>.
 */
int VECTOR_PTR_GetSize(VECTOR_PTR* vec);

/*! FUNCTION:  VECTOR_PTR_SetSize()
 *  SYNOPSIS:  Set utilized length of <vec>
 *  RETURN:    Pointer to location to <vec> idx.
 */
STATUS_FLAG
VECTOR_PTR_SetSize(VECTOR_PTR* vec, size_t size);

/*! FUNCTION:  VECTOR_PTR_GetSizeAlloc()
 *  SYNOPSIS:  Get allocated length of <vec>.
 */
int VECTOR_PTR_GetSizeAlloc(VECTOR_PTR* vec);

/*! FUNCTION:  VECTOR_PTR_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array.
 *             Returns index of the first occurance found.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of occurance of <val>.
 *             Return -1 if <val> is not found.
 */
int VECTOR_PTR_Search(VECTOR_PTR* vec, PTR val);

/*! FUNCTION:  VECTOR_PTR_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array.
 *             Returns index of the first occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of first occurance of <val>.
 *             Return -1 if <val> is not found.
 */
int VECTOR_PTR_Search_First(VECTOR_PTR* vec, PTR val);

/*! FUNCTION:  VECTOR_PTR_Search_Last()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array.
 *             Returns index of the last occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of last occurance of <val>.
 *             Return -1 if <val> is not found.
 */
int VECTOR_PTR_Search_Last(VECTOR_PTR* vec, PTR val);

/*! FUNCTION:  VECTOR_PTR_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality,
 *             pos if <vec_A> > <vec_B>,
 *             neg if <vec_A> < <vec_B>.
 */
int VECTOR_PTR_Compare(VECTOR_PTR* vec_A, VECTOR_PTR* vec_B);

/*! FUNCTION:  VECTOR_PTR_Sort()
 *  SYNOPSIS:  Sorts <vec> data array in ascending order. In-place.
 */
STATUS_FLAG
VECTOR_PTR_Sort(VECTOR_PTR* vec);

/*! FUNCTION:  VECTOR_PTR_Sort_Sub()
 *  SYNOPSIS:  Subcall - Sorts <vec> data array in ascending order on range
 * (beg,end]. Uses quicksort until length of subarray falls below threshold,
 * then selection sort.
 */
STATUS_FLAG
VECTOR_PTR_Sort_Sub(VECTOR_PTR* vec, int beg, int end);

/*! FUNCTION:  VECTOR_PTR_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs selection sort on <vec> data array in ascending
 * order on range (beg,end].
 */
STATUS_FLAG
VECTOR_PTR_Sort_Sub_Selectsort(VECTOR_PTR* vec, int beg, int end);

/*! FUNCTION:  VECTOR_PTR_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Subcall - Runs quicksort on <vec> data array in ascending order
 * on range (beg,end].
 */
STATUS_FLAG
VECTOR_PTR_Sort_Sub_Quicksort(VECTOR_PTR* vec, int beg, int end);

/*! FUNCTION:  VECTOR_PTR_Op()
 *  SYNOPSIS:  Perform element-wise unary operation <op>(PTR data) to each cell
 * in <vec_in> and puts it in <vec_out>. Returns a pointer to <vec_out>.
 *             <vec_out> will be resized to size of <vec_in>.
 *             <vec_in> and <vec_out> can be the same vector. If <vec_out> is
 * NULL, new vector will be created.
 */
VECTOR_PTR* VECTOR_PTR_Op(
    VECTOR_PTR* vec_out,  /* input vector */
    VECTOR_PTR* vec_in,   /* output vector (can be input vector) */
    PTR (*op)(PTR data)); /* unary operation function */

/*! FUNCTION:  VECTOR_PTR_Op()
 *  SYNOPSIS:  Perform element-wise binary operation <op>(PTR data_1, PTR
 * data_2) to each cell in <vec_in_1, vec_in_2> and puts it in <vec_out>.
 *             <vec_in_1> and <vec_in_2> must be the same size.
 *             Returns a pointer to <vec_out>.
 *             <vec_out> will be resized to size of <vec_in>.
 *             <vec_in> and <vec_out> can be the same vector. If <vec_out> is
 * NULL, new vector will be created.
 */
VECTOR_PTR* VECTOR_PTR_BinOp(
    VECTOR_PTR* vec_out,                /* output vector */
    VECTOR_PTR* vec_in_1,               /* first input vector */
    VECTOR_PTR* vec_in_2,               /* second input vector */
    PTR (*op)(PTR data_1, PTR data_2)); /* binary operation function */

/*! FUNCTION:  VECTOR_PTR_Swap()
 *  SYNOPSIS:  Swaps the values of <vec> at indexes <i> and <j>.
 */
STATUS_FLAG
VECTOR_PTR_Swap(VECTOR_PTR* vec, int i, int j);

/*! FUNCTION:  VECTOR_PTR_Reverse()
 *  SYNOPSIS:  Reverse the ordering of array.
 */
STATUS_FLAG
VECTOR_PTR_Reverse(VECTOR_PTR* vec);

/*! FUNCTION:  VECTOR_PTR_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
STATUS_FLAG
VECTOR_PTR_Dump(VECTOR_PTR* vec, FILE* fp);

/*! FUNCTION:  VECTOR_PTR_Dump_byOpt()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer. Non-optimized.
 */
STATUS_FLAG
VECTOR_PTR_Dump_byOpt(VECTOR_PTR* vec, STR delim, STR header, FILE* fp);

/*! FUNCTION:  VECTOR_PTR_Unit_Test()
 *  SYNOPSIS:  Perform unit test for VECTOR_PTR.
 */
STATUS_FLAG
VECTOR_PTR_Unit_Test();

#endif
