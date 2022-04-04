/*******************************************************************************
 *  - FILE:   vector_trace.c
 *  - DESC:    VECTOR_TRACE Object Functions.
 *             Provides template for building vector classes.
 *             Run "scripts/builder-helper/build_vector_classes_from_template" to update.
 *             Requires data primitive to have TRACE_Compare().
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>

/* local imports */
#include "../../objects/structs.h"
#include "../../utilities/_utilities.h"
#include "../../objects/_objects.h"

/* header */
#include "_vectors.h"
#include "vector_trace.h"

/*! FUNCTION:  VECTOR_TRACE_Create()
 *  SYNOPSIS:  Create new VECTOR_TRACE object and returns pointer.
 */
VECTOR_TRACE*
VECTOR_TRACE_Create() {
  const int init_size = VECTOR_INIT_SIZE;
  return VECTOR_TRACE_Create_by_Size(init_size);
}

/*! FUNCTION:  VECTOR_TRACE_Create()
 *  SYNOPSIS:  Create new VECTOR_TRACE object at specific size and returns pointer.
 */
VECTOR_TRACE*
VECTOR_TRACE_Create_by_Size(size_t size) {
  VECTOR_TRACE* vec;
  vec = ERROR_malloc(sizeof(VECTOR_TRACE));

  vec->data = NULL;
  vec->N = 0;
  vec->Nalloc = 0;

  VECTOR_TRACE_Resize(vec, size);

  return vec;
}

/*! FUNCTION:  VECTOR_TRACE_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_TRACE.
 */
VECTOR_TRACE*
VECTOR_TRACE_Destroy(VECTOR_TRACE* vec) {
  if (vec == NULL)
    return NULL;

  /* NOTE: This loop prevents memory leaks with data types which allocate dynamic memory
   *       Hopefully, this can be safely used in template, as this loop should be
   *       inlined and optimized out by non-dynamic data types(?)
   */
  // for (int i = 0; i < vec->N; i++) {
  //    VEC_X( vec, i ) = TRACE_Destroy( VEC_X( vec, i ) );
  // }

  vec->data = ERROR_free(vec->data);
  vec = ERROR_free(vec);

  return NULL;
}

/*! FUNCTION:  VECTOR_TRACE_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_TRACE object by resetting size counter (no realloc) .
 */
STATUS_FLAG
VECTOR_TRACE_Reuse(VECTOR_TRACE* vec) {
  /* NOTE: This loop prevents memory leaks with data types which allocate dynamic memory
   *       Hopefully, this can be safely used in template, as this loop should be
   *       inlined and optimized out by non-dynamic data types(?)
   */
  // for (int i = 0; i < vec->N; i++) {
  //    VEC_X( vec, i ) = TRACE_Destroy( VEC_X( vec, i ) );
  // }

  vec->N = 0;
}

/* TODO: Can I do this without dynamic allocation? */
/*! FUNCTION:  VECTOR_TRACE_WrapArray().
 *  SYNOPSIS:  Wraps <array> in a <vector> struct, allocates structures and returns pointer.
 *             <length> is the entire array size, <occupied> is the amount of data in use
 *             WARNING: Limited operation support. Do not use:
 *                      _SetSize()     *Unless size is less than total array size
 *                      _Pushback()    *Use _Push()
 *                      _Popback()     *Use _Pop()
 */
VECTOR_TRACE*
VECTOR_TRACE_WrapArray(TRACE* array,
                       size_t length,
                       size_t occupied) {
  VECTOR_TRACE* vec;
  vec = ERROR_malloc(sizeof(VECTOR_TRACE));

  vec->data = array;
  vec->N = occupied;
  vec->Nalloc = length;

  return vec;
}

/*! FUNCTION:  VECTOR_TRACE_UnwrapArray()
 *  SYNOPSIS:  Unwraps <array> from <vector> and return it.  Frees all <vector> related data.
 */
TRACE*
VECTOR_TRACE_UnwrapArray(VECTOR_TRACE* vec,
                         size_t size) {
  TRACE* array = vec->data;
  vec = ERROR_malloc(sizeof(VECTOR_TRACE));

  return array;
}

/*! FUNCTION:  VECTOR_TRACE_GetArray()
 *  SYNOPSIS:  Get <data> array from <vec>.
 */
inline TRACE*
VECTOR_TRACE_GetArray(VECTOR_TRACE* vec) {
  return vec->data;
}

/*! FUNCTION:  VECTOR_TRACE_Fill()
 *  SYNOPSIS:  Fill VECTOR_TRACE object with val.
 */
STATUS_FLAG
VECTOR_TRACE_Fill(VECTOR_TRACE* vec,
                  TRACE val) {
  for (int i = 0; i < vec->N; i++) {
    VEC_X(vec, i) = TRACE_Create(val);
  }
}

/*! FUNCTION:  VECTOR_TRACE_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object.
 *             Creates new VECTOR_TRACE for <dest> if <dest> is NULL.
 */
VECTOR_TRACE*
VECTOR_TRACE_Copy(VECTOR_TRACE* dest,
                  VECTOR_TRACE* src) {
  if (dest == NULL) {
    dest = VECTOR_TRACE_Create();
  }
  /* allocate variable-sized data */
  VECTOR_TRACE_GrowTo(dest, src->N);
  /* copy variable-sized data */
  for (int i = 0; i < src->N; i++) {
    *VECTOR_TRACE_GetX(dest, i) = TRACE_Create(VEC_X(src, i));
  }
  /* copy base data */
  dest->N = src->N;

  return dest;
}

/*! FUNCTION:  VECTOR_TRACE_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>.
 */
STATUS_FLAG
VECTOR_TRACE_Resize(VECTOR_TRACE* vec,
                    size_t size) {
  vec->data = ERROR_realloc(vec->data, sizeof(TRACE) * size);
  vec->Nalloc = size;
}

/*! FUNCTION:  VECTOR_TRACE_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>.
 */
STATUS_FLAG
VECTOR_TRACE_GrowTo(VECTOR_TRACE* vec,
                    size_t size) {
  if (vec->Nalloc < size) {
    VECTOR_TRACE_Resize(vec, size);
  }
}

/*! FUNCTION:  VECTOR_TRACE_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return value of data.
 *  RETURN:    Return data at <idx>.
 */
inline TRACE
VECTOR_TRACE_Get(VECTOR_TRACE* vec,
                 int idx) {
/* if debugging, do edgebound checks */
#if SAFE
  if (idx >= vec->N || idx < 0) {
    fprintf(stderr, "ERROR: VECTOR_TRACE access out-of-bounds.\n");
    fprintf(stderr, "dim: (%ld/%ld), access: %d\n", vec->N, vec->Nalloc, idx);
    ERRORCHECK_exit(EXIT_FAILURE);
  }
#endif

  return (vec->data[idx]);
}

/*! FUNCTION:  VECTOR_TRACE_GetX()
 *  SYNOPSIS:  Get reference to data from <vec> at the <idx> position in array, and return pointer to data.
 *             Warning: Out-of-Bounds only checked in DEBUG and SAFE.
 *  RETURN:    Pointer to location to <vec> idx.
 */
inline TRACE*
VECTOR_TRACE_GetX(VECTOR_TRACE* vec,
                  int idx) {
/* if debugging, do edgebound checks */
#if SAFE
  if (idx >= vec->N || idx < 0) {
    fprintf(stderr, "ERROR: VECTOR_TRACE access out-of-bounds.\n");
    fprintf(stderr, "dim: (%ld/%ld), access: %d\n", vec->N, vec->Nalloc, idx);
    ERRORCHECK_exit(EXIT_FAILURE);
  }
#endif

  return &(vec->data[idx]);
}

/*! FUNCTION:  VECTOR_TRACE_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
STATUS_FLAG
VECTOR_TRACE_Set(VECTOR_TRACE* vec,
                 int idx,
                 TRACE val) {
/* if debugging, do edgebound checks */
#if SAFE
  if (idx >= vec->N || idx < 0) {
    fprintf(stderr, "ERROR: VECTOR_TRACE access out-of-bounds.\n");
    fprintf(stderr, "dim: (%ld/%ld), access: %d\n", vec->N, vec->Nalloc, idx);
    ERRORCHECK_exit(EXIT_FAILURE);
  }
#endif

  VEC_X(vec, idx) = TRACE_Destroy(VEC_X(vec, idx));
  VEC_X(vec, idx) = TRACE_Create(val);

  return STATUS_SUCCESS;
}

/*! FUNCTION:  VECTOR_TRACE_Insert()
 *  SYNOPSIS:  Overwrite data from <vec> at the <idx> position in array to <val>. Deletes present value.
 */
STATUS_FLAG
VECTOR_TRACE_Insert(VECTOR_TRACE* vec,
                    int idx,
                    TRACE val) {
/* if debugging, do edgebound checks */
#if SAFE
  if (idx >= vec->N || idx < 0) {
    fprintf(stderr, "ERROR: VECTOR_TRACE access out-of-bounds.\n");
    fprintf(stderr, "dim: (%ld/%ld), access: %d\n", vec->N, vec->Nalloc, idx);
    ERRORCHECK_exit(EXIT_FAILURE);
  }
#endif

  VEC_X(vec, idx) = TRACE_Destroy(VEC_X(vec, idx));
  VEC_X(vec, idx) = TRACE_Create(val);

  return STATUS_SUCCESS;
}

/*! FUNCTION:  VECTOR_TRACE_Delete()
 *  SYNOPSIS:  Overwrite data from <vec> at the <idx> position in array to <val>. Deletes present value.
 */
STATUS_FLAG
VECTOR_TRACE_Delete(VECTOR_TRACE* vec,
                    int idx) {
/* if debugging, do edgebound checks */
#if SAFE
  if (idx >= vec->N || idx < 0) {
    fprintf(stderr, "ERROR: VECTOR_TRACE access out-of-bounds.\n");
    fprintf(stderr, "dim: (%ld/%ld), access: %d\n", vec->N, vec->Nalloc, idx);
    ERRORCHECK_exit(EXIT_FAILURE);
  }
#endif

  int N = VECTOR_TRACE_GetSize(vec);
  VEC_X(vec, idx) = TRACE_Destroy(VEC_X(vec, idx));
  VEC_X(vec, idx) = VEC_X(vec, N - 1);
  // VEC_X( vec, N - 1 )   = TRACE_Empty();
  vec->N -= 1;

  return STATUS_SUCCESS;
}

/*! FUNCTION:  VECTOR_TRACE_Push()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array.
 *             WARNING: Does not handle resizing or check for out-of-bounds.
 *                      For that, use Pushback().
 */
inline STATUS_FLAG
VECTOR_TRACE_Push(VECTOR_TRACE* vec,
                  TRACE val) {
  /* NOTE: This push() creates another copy of the data to store in vector (in the case of dynamically allocated data) */
  VEC_X(vec, vec->N) = TRACE_Create(val);
  vec->N++;
}

/*! FUNCTION:  VECTOR_TRACE_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             Resize array if array is full.
 */
inline STATUS_FLAG
VECTOR_TRACE_Pushback(VECTOR_TRACE* vec,
                      TRACE val) {
  VECTOR_TRACE_Push(vec, val);

  /* if array is full, resize */
  if (vec->N >= vec->Nalloc - 1) {
    VECTOR_TRACE_Resize(vec, vec->N * 2);
  }
}

/*! FUNCTION:  VECTOR_TRACE_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, remove data, and return data.
 */
inline TRACE
VECTOR_TRACE_Pop(VECTOR_TRACE* vec) {
  TRACE data = VECTOR_TRACE_Get(vec, vec->N - 1);
  vec->N -= 1;

  return data;
}

/*! FUNCTION:  VECTOR_TRACE_Popback()
 *  SYNOPSIS:  Pop data from the end of <vec> data array and return data.
 *             Resize if array is less than half full.
 */
inline TRACE
VECTOR_TRACE_Popback(VECTOR_TRACE* vec) {
  TRACE data = VECTOR_TRACE_Pop(vec);

  /* if array is less than half used, resize */
  if (vec->N < vec->Nalloc / 2) {
    VECTOR_TRACE_Resize(vec, vec->N / 2);
  }

  return data;
}

/*! FUNCTION:  VECTOR_TRACE_Append()
 *  SYNOPSIS:  Push <append> data array of length <L> onto the end of <vec> data array.
 */
inline STATUS_FLAG
VECTOR_TRACE_Append(VECTOR_TRACE* vec,
                    TRACE* append,
                    size_t L) {
  size_t N_new;
  N_new = vec->N + L;

  /* resize array */
  if (vec->Nalloc < N_new) {
    VECTOR_TRACE_Resize(vec, N_new);
  }

  /* copy data over */
  for (int i = 0; i < L; i++) {
    VECTOR_TRACE_Push(vec, append[i]);
  }
}

/*! FUNCTION:  VECTOR_TRACE_GetSize()
 *  SYNOPSIS:  Get utilized length of <vec>.
 */
inline int
VECTOR_TRACE_GetSize(VECTOR_TRACE* vec) {
  return vec->N;
}

/*! FUNCTION:  VECTOR_TRACE_SetSize()
 *  SYNOPSIS:  Set utilized length of <vec>.
 *             Will allocate memory if necessary.
 */
inline STATUS_FLAG
VECTOR_TRACE_SetSize(VECTOR_TRACE* vec,
                     size_t size) {
  VECTOR_TRACE_GrowTo(vec, size);
  vec->N = size;
}

/*! FUNCTION:  VECTOR_TRACE_GetSizeAlloc()
 *  SYNOPSIS:  Get allocated length of <vec>.
 */
inline int
VECTOR_TRACE_GetSizeAlloc(VECTOR_TRACE* vec) {
  return vec->Nalloc;
}

/*! FUNCTION:  VECTOR_TRACE_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array.
 *             Returns index of the first occurance found.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of occurance of <val>.
 *             Return -1 if <val> is not found.
 */
int VECTOR_TRACE_Search(VECTOR_TRACE* vec,
                        TRACE val) {
  int N = VECTOR_TRACE_GetSize(vec);
  int idx = (N / 2);
  int found = -1;

  for (int i = N / 4; i >= 1; i /= 2) {
    int cmp = TRACE_Compare(val, vec->data[i]);

    if (cmp > 0) {
      idx += i;
    } else if (cmp < 0) {
      idx -= i;
    } else {
      found = idx;
      return found;
    }
  }
  return found;
}

/*! FUNCTION:  VECTOR_TRACE_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array.
 *             Returns index of the first occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of first occurance of <val>.
 *             Return -1 if <val> is not found.
 */
int VECTOR_TRACE_Search_First(VECTOR_TRACE* vec,
                              TRACE val) {
  int N = VECTOR_TRACE_GetSize(vec);
  int idx = (N / 2);
  int found = -1;

  for (int i = N / 4; i >= 1; i /= 2) {
    int cmp = TRACE_Compare(val, vec->data[i]);

    if (cmp > 0) {
      idx += i;
    } else if (cmp < 0) {
      idx -= i;
    } else {
      found = idx;
      idx -= i;
    }
  }
  return found;
}

/*! FUNCTION:  VECTOR_TRACE_Search_Last()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array.
 *             Returns index of the last occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of last occurance of <val>.
 *             Return -1 if <val> is not found.
 */
int VECTOR_TRACE_Search_Last(VECTOR_TRACE* vec,
                             TRACE val) {
  int N = VECTOR_TRACE_GetSize(vec);
  int idx = (N / 2);
  int found = -1;

  for (int i = N / 4; i >= 1; i /= 2) {
    int cmp = TRACE_Compare(val, vec->data[i]);

    if (cmp > 0) {
      idx += i;
    } else if (cmp < 0) {
      idx -= i;
    } else {
      found = idx;
      idx += i;
    }
  }
  return found;
}

/*! FUNCTION:  VECTOR_TRACE_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality,
 *             pos if <vec_A> > <vec_B>,
 *             neg if <vec_A> < <vec_B>.
 */
int VECTOR_TRACE_Compare(VECTOR_TRACE* vec_A,
                         VECTOR_TRACE* vec_B) {
  for (int i = 0; i < vec_A->N; i++) {
    if (TRACE_Compare(vec_A->data[i], vec_B->data[i]) != 0) {
      return TRACE_Compare(vec_A->data[i], vec_B->data[i]);
    }
  }
  return 0;
}

/*! FUNCTION:  VECTOR_TRACE_Sort()
 *  SYNOPSIS:  Sorts <vec> data array in ascending order. In-place.
 */
STATUS_FLAG
VECTOR_TRACE_Sort(VECTOR_TRACE* vec) {
  int N = VECTOR_TRACE_GetSize(vec);
  qsort(vec->data, N, sizeof(TRACE), TRACE_CompareTo);

  // VECTOR_TRACE_Sort_Sub( vec, 0, N );

#if DEBUG
  {
    for (int i = 0; i < N - 1; i++) {
      TRACE cur = vec->data[i];
      TRACE nxt = vec->data[i + 1];
      char s_cur[50];
      char s_nxt[50];
      int cmp = TRACE_Compare(cur, nxt);
      if ((cmp <= 0) == false) {
        fprintf(stderr, "ERROR: bad sort. %d, %d v %d: %s vs %s\n",
                cmp, i, i + 1, TRACE_ToString(cur, s_cur), TRACE_ToString(nxt, s_nxt));
        ERRORCHECK_exit(EXIT_FAILURE);
      }
    }
  }
#endif
}

/*! FUNCTION:  VECTOR_TRACE_Sort_Sub()
 *  SYNOPSIS:  Sorts subarray of <vec> data in ascending order on range (beg,end].
 *             Uses quicksort until length of subarray falls below threshold, then selection sort.
 */
STATUS_FLAG
VECTOR_TRACE_Sort_Sub(VECTOR_TRACE* vec,
                      int beg,
                      int end) {
  const int begin_select_sort = INT_MAX;
  int N = VECTOR_TRACE_GetSize(vec);
  if (N <= 1) {
    return STATUS_SUCCESS;
  }

  int size = end - beg;
  /* run selection sort if below threshold */
  if (size <= begin_select_sort) {
    VECTOR_TRACE_Sort_Sub_Selectsort(vec, beg, end);
  }
  /* otherwise run quiksort */
  else {
    VECTOR_TRACE_Sort_Sub_Quicksort(vec, beg, end);
  }

  return STATUS_SUCCESS;
}

/*! FUNCTION:  VECTOR_TRACE_Sort_Sub_Selectsort()
 *  SYNOPSIS:  Selection Sorts subarray of <vec> data in ascending order on range (beg,end].
 */
STATUS_FLAG
VECTOR_TRACE_Sort_Sub_Selectsort(VECTOR_TRACE* vec,
                                 int beg,
                                 int end) {
  /* find the minimum element of remaining unsorted list */
  for (int i = beg; i < end; i++) {
    /* initial minimum value found */
    int min_idx = i;
    TRACE min_val = vec->data[i];
    for (int j = i + 1; j < end; j++) {
      /* if new minimum found, update value and index */
      int cmp = TRACE_Compare(min_val, vec->data[j]);
      if (cmp > 0) {
        min_idx = j;
        min_val = vec->data[j];
      }
    }
    /* swap new minimum to left-most position */
    VECTOR_TRACE_Swap(vec, i, min_idx);
  }
}

/*! FUNCTION:  VECTOR_TRACE_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Quick Sorts subarray of <vec> data in ascending order on range (beg,end].
 */
STATUS_FLAG
VECTOR_TRACE_Sort_Sub_Quicksort(VECTOR_TRACE* vec,
                                int beg,
                                int end) {
  /* partition pointers */
  int r_idx = beg + 1;
  int l_idx = end - 1;
  TRACE* rhs = &(vec->data[beg + 1]);
  TRACE* lhs = &(vec->data[end - 1]);

  /* select random pivot value */
  int range = end - beg;
  int pivot_idx = RNG_INT_Range(beg, end);
  TRACE pivot_val = vec->data[pivot_idx];
  VECTOR_TRACE_Swap(vec, pivot_idx, beg);

  /* partition on pivot */
  while (l_idx <= r_idx) {
    /* find next right partition element that is less than pivot element */
    while ((l_idx <= r_idx) && (TRACE_Compare(pivot_val, vec->data[r_idx]) < 0)) {
      r_idx--;
    }
    /* find next left partition element that is greater than pivot element */
    while ((l_idx <= r_idx) && (TRACE_Compare(pivot_val, vec->data[l_idx]) >= 0)) {
      l_idx++;
    }
    /* if left and right index have not crossed, then swap elements */
    if (l_idx <= r_idx) {
      VECTOR_TRACE_Swap(vec, l_idx, r_idx);
    }
  }
  /* move partition element to barrier between left and right index */
  VECTOR_TRACE_Swap(vec, beg, r_idx);
  /* sort both partitions (omit partition element) */
  VECTOR_TRACE_Sort_Sub(vec, beg, r_idx);
  VECTOR_TRACE_Sort_Sub(vec, r_idx + 1, end);
}

/*! FUNCTION:  VECTOR_TRACE_Op()
 *  SYNOPSIS:  Perform element-wise unary operation <op> to each cell in <vec_in> and puts it in <vec_out>.
 *             Returns a pointer to <vec_out>.
 *             <vec_out> will be resized to size of <vec_in>.
 *             <vec_in> and <vec_out> can be the same vector. If <vec_out> is NULL, new vector will be created.
 */
inline VECTOR_TRACE*
VECTOR_TRACE_Op(VECTOR_TRACE* vec_out,   /* input vector */
                VECTOR_TRACE* vec_in,    /* output vector (can be input vector) */
                TRACE (*op)(TRACE data)) /* unary operation */
{
  int N = VECTOR_TRACE_GetSize(vec_in);
  if (vec_out == NULL) {
    VECTOR_TRACE_Create_by_Size(N);
  }
  VECTOR_TRACE_SetSize(vec_out, N);

  for (int i = 0; i < N; i++) {
    *VECTOR_TRACE_GetX(vec_out, i) = op(*VECTOR_TRACE_GetX(vec_in, i));
  }

  return vec_out;
}

/*! FUNCTION:  VECTOR_TRACE_Op()
 *  SYNOPSIS:  Perform element-wise binary operation <op>(TRACE data_1, TRACE data_2) to each cell in <vec_in_1, vec_in_2> and puts it in <vec_out>.
 *             <vec_in_1> and <vec_in_2> must be the same size.
 *             Returns a pointer to <vec_out>.
 *             <vec_out> will be resized to size of <vec_in>.
 *             <vec_in> and <vec_out> can be the same vector. If <vec_out> is NULL, new vector will be created.
 */
inline VECTOR_TRACE*
VECTOR_TRACE_BinOp(VECTOR_TRACE* vec_out,                   /* output vector */
                   VECTOR_TRACE* vec_in_1,                  /* first input vector */
                   VECTOR_TRACE* vec_in_2,                  /* second input vector */
                   TRACE (*op)(TRACE data_1, TRACE data_2)) /* binary operation */
{
  int N = VECTOR_TRACE_GetSize(vec_in_1);
  int N_2 = VECTOR_TRACE_GetSize(vec_in_2);
  if (N != N_2) {
    fprintf(stderr, "ERROR: <vec_in_1> and <vec_in_2> are not the same dimensions.\n");
    ERRORCHECK_exit(EXIT_FAILURE);
  }
  if (vec_out == NULL) {
    VECTOR_TRACE_Create_by_Size(N);
  }
  VECTOR_TRACE_SetSize(vec_out, N);

  for (int i = 0; i < N; i++) {
    *VECTOR_TRACE_GetX(vec_out, i) = op(*VECTOR_TRACE_GetX(vec_in_1, i), *VECTOR_TRACE_GetX(vec_in_2, i));
  }

  return vec_out;
}

/*! FUNCTION:  VECTOR_TRACE_Swap()
 *  SYNOPSIS:  Swaps the values of <vec> at indexes <i> and <j>.
 *             Warning: Only checks for Out-of-Bounds when in DEBUG.
 */
inline STATUS_FLAG
VECTOR_TRACE_Swap(VECTOR_TRACE* vec,
                  int i,
                  int j) {
  TRACE swap = VECTOR_TRACE_Get(vec, i);
  vec->data[i] = vec->data[j];
  vec->data[j] = swap;
}

/*! FUNCTION:  VECTOR_TRACE_Reverse()
 *  SYNOPSIS:  Reverse the ordering of array.
 */
inline STATUS_FLAG VECTOR_TRACE_Reverse(VECTOR_TRACE* vec) {
  int N = VECTOR_TRACE_GetSize(vec);

  for (int i = 0; i < (N / 2); i++) {
    VECTOR_TRACE_Swap(vec, i, (N - 1) - i);
  }
}

/*! FUNCTION:  VECTOR_TRACE_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer. Non-optimized.
 */
STATUS_FLAG
VECTOR_TRACE_Dump(VECTOR_TRACE* vec,
                  FILE* fp) {
  VECTOR_TRACE_Dump_byOpt(vec, "\n", "VECTOR BOUND", fp);
}

/*! FUNCTION:  VECTOR_TRACE_Dump_byOpt()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
STATUS_FLAG
VECTOR_TRACE_Dump_byOpt(VECTOR_TRACE* vec,
                        STR delim,
                        STR header,
                        FILE* fp) {
  /* stringification of template object */
  char s[50];
  const char* pad = " ";

  fprintf(fp, "%s: ", header);
  fprintf(fp, "[ ");
  for (int i = 0; i < vec->N; i++) {
    fprintf(fp, "%s%s%s", TRACE_ToString(vec->data[i], s), delim, pad);
  }
  if (vec->N >= 1) {
    fprintf(fp, "%s%s", TRACE_ToString(vec->data[vec->N - 1], s), pad);
  }
  fprintf(fp, "]\n");
}

/*! FUNCTION:  VECTOR_TRACE_UnitTest()
 *  SYNOPSIS:  Perform unit test for VECTOR_TRACE.
 */
STATUS_FLAG
VECTOR_TRACE_UnitTest() {
  VECTOR_TRACE* vec = VECTOR_TRACE_Create();

  VECTOR_TRACE_Destroy(vec);
}
