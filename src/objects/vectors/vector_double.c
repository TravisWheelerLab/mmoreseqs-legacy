/*******************************************************************************
 *  - FILE:   vector_double.c
 *  - DESC:    VECTOR_DBL Object Functions.
 *             Provides template for building vector classes.
 *             Run "scripts/builder-helper/build_vector_classes_from_template" to update.
 *             Requires data primitive to have DBL_Compare().
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
#include "vector_double.h"

/*! FUNCTION:  VECTOR_DBL_Create()
 *  SYNOPSIS:  Create new VECTOR_DBL object and returns pointer.
 */
VECTOR_DBL*
VECTOR_DBL_Create() {
  const int init_size = VECTOR_INIT_SIZE;
  return VECTOR_DBL_Create_by_Size(init_size);
}

/*! FUNCTION:  VECTOR_DBL_Create()
 *  SYNOPSIS:  Create new VECTOR_DBL object at specific size and returns pointer.
 */
VECTOR_DBL*
VECTOR_DBL_Create_by_Size(size_t size) {
  VECTOR_DBL* vec;
  vec = ERROR_malloc(sizeof(VECTOR_DBL));

  vec->data = NULL;
  vec->N = 0;
  vec->Nalloc = 0;

  VECTOR_DBL_Resize(vec, size);

  return vec;
}

/*! FUNCTION:  VECTOR_DBL_Destroy()
 *  SYNOPSIS:  Frees all data associated with VECTOR_DBL.
 */
VECTOR_DBL*
VECTOR_DBL_Destroy(VECTOR_DBL* vec) {
  if (vec == NULL)
    return NULL;

  /* NOTE: This loop prevents memory leaks with data types which allocate dynamic memory
   *       Hopefully, this can be safely used in template, as this loop should be
   *       inlined and optimized out by non-dynamic data types(?)
   */
  // for (int i = 0; i < vec->N; i++) {
  //    VEC_X( vec, i ) = DBL_Destroy( VEC_X( vec, i ) );
  // }

  vec->data = ERROR_free(vec->data);
  vec = ERROR_free(vec);

  return NULL;
}

/*! FUNCTION:  VECTOR_DBL_Reuse()
 *  SYNOPSIS:  Reuse VECTOR_DBL object by resetting size counter (no realloc) .
 */
STATUS_FLAG
VECTOR_DBL_Reuse(VECTOR_DBL* vec) {
  /* NOTE: This loop prevents memory leaks with data types which allocate dynamic memory
   *       Hopefully, this can be safely used in template, as this loop should be
   *       inlined and optimized out by non-dynamic data types(?)
   */
  // for (int i = 0; i < vec->N; i++) {
  //    VEC_X( vec, i ) = DBL_Destroy( VEC_X( vec, i ) );
  // }

  vec->N = 0;
}

/* TODO: Can I do this without dynamic allocation? */
/*! FUNCTION:  VECTOR_DBL_WrapArray().
 *  SYNOPSIS:  Wraps <array> in a <vector> struct, allocates structures and returns pointer.
 *             <length> is the entire array size, <occupied> is the amount of data in use
 *             WARNING: Limited operation support. Do not use:
 *                      _SetSize()     *Unless size is less than total array size
 *                      _Pushback()    *Use _Push()
 *                      _Popback()     *Use _Pop()
 */
VECTOR_DBL*
VECTOR_DBL_WrapArray(DBL* array,
                     size_t length,
                     size_t occupied) {
  VECTOR_DBL* vec;
  vec = ERROR_malloc(sizeof(VECTOR_DBL));

  vec->data = array;
  vec->N = occupied;
  vec->Nalloc = length;

  return vec;
}

/*! FUNCTION:  VECTOR_DBL_UnwrapArray()
 *  SYNOPSIS:  Unwraps <array> from <vector> and return it.  Frees all <vector> related data.
 */
DBL* VECTOR_DBL_UnwrapArray(VECTOR_DBL* vec,
                            size_t size) {
  DBL* array = vec->data;
  vec = ERROR_malloc(sizeof(VECTOR_DBL));

  return array;
}

/*! FUNCTION:  VECTOR_DBL_GetArray()
 *  SYNOPSIS:  Get <data> array from <vec>.
 */
inline DBL*
VECTOR_DBL_GetArray(VECTOR_DBL* vec) {
  return vec->data;
}

/*! FUNCTION:  VECTOR_DBL_Fill()
 *  SYNOPSIS:  Fill VECTOR_DBL object with val.
 */
STATUS_FLAG
VECTOR_DBL_Fill(VECTOR_DBL* vec,
                DBL val) {
  for (int i = 0; i < vec->N; i++) {
    VEC_X(vec, i) = DBL_Create(val);
  }
}

/*! FUNCTION:  VECTOR_DBL_Copy()
 *  SYNOPSIS:  Create deep copy of <src> object.
 *             Creates new VECTOR_DBL for <dest> if <dest> is NULL.
 */
VECTOR_DBL*
VECTOR_DBL_Copy(VECTOR_DBL* dest,
                VECTOR_DBL* src) {
  if (dest == NULL) {
    dest = VECTOR_DBL_Create();
  }
  /* allocate variable-sized data */
  VECTOR_DBL_GrowTo(dest, src->N);
  /* copy variable-sized data */
  for (int i = 0; i < src->N; i++) {
    *VECTOR_DBL_GetX(dest, i) = DBL_Create(VEC_X(src, i));
  }
  /* copy base data */
  dest->N = src->N;

  return dest;
}

/*! FUNCTION:  VECTOR_DBL_Resize()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>.
 */
STATUS_FLAG
VECTOR_DBL_Resize(VECTOR_DBL* vec,
                  size_t size) {
  vec->data = ERROR_realloc(vec->data, sizeof(DBL) * size);
  vec->Nalloc = size;
}

/*! FUNCTION:  VECTOR_DBL_GrowTo()
 *  SYNOPSIS:  Reallocate <vec> data array to length of <size>,
 *             only if current array length is less than <size>.
 */
STATUS_FLAG
VECTOR_DBL_GrowTo(VECTOR_DBL* vec,
                  size_t size) {
  if (vec->Nalloc < size) {
    VECTOR_DBL_Resize(vec, size);
  }
}

/*! FUNCTION:  VECTOR_DBL_Get()
 *  SYNOPSIS:  Get data from <vec> at the <idx> position in array, and return value of data.
 *  RETURN:    Return data at <idx>.
 */
inline DBL
VECTOR_DBL_Get(VECTOR_DBL* vec,
               int idx) {
/* if debugging, do edgebound checks */
#if SAFE
  if (idx >= vec->N || idx < 0) {
    fprintf(stderr, "ERROR: VECTOR_DBL access out-of-bounds.\n");
    fprintf(stderr, "dim: (%ld/%ld), access: %d\n", vec->N, vec->Nalloc, idx);
    ERRORCHECK_exit(EXIT_FAILURE);
  }
#endif

  return (vec->data[idx]);
}

/*! FUNCTION:  VECTOR_DBL_GetX()
 *  SYNOPSIS:  Get reference to data from <vec> at the <idx> position in array, and return pointer to data.
 *             Warning: Out-of-Bounds only checked in DEBUG and SAFE.
 *  RETURN:    Pointer to location to <vec> idx.
 */
inline DBL*
VECTOR_DBL_GetX(VECTOR_DBL* vec,
                int idx) {
/* if debugging, do edgebound checks */
#if SAFE
  if (idx >= vec->N || idx < 0) {
    fprintf(stderr, "ERROR: VECTOR_DBL access out-of-bounds.\n");
    fprintf(stderr, "dim: (%ld/%ld), access: %d\n", vec->N, vec->Nalloc, idx);
    ERRORCHECK_exit(EXIT_FAILURE);
  }
#endif

  return &(vec->data[idx]);
}

/*! FUNCTION:  VECTOR_DBL_Set()
 *  SYNOPSIS:  Set data from <vec> at the <idx> position in array to <val>.
 *             Warning: Out-of-Bounds only checked in DEBUG.
 */
STATUS_FLAG
VECTOR_DBL_Set(VECTOR_DBL* vec,
               int idx,
               DBL val) {
/* if debugging, do edgebound checks */
#if SAFE
  if (idx >= vec->N || idx < 0) {
    fprintf(stderr, "ERROR: VECTOR_DBL access out-of-bounds.\n");
    fprintf(stderr, "dim: (%ld/%ld), access: %d\n", vec->N, vec->Nalloc, idx);
    ERRORCHECK_exit(EXIT_FAILURE);
  }
#endif

  VEC_X(vec, idx) = DBL_Destroy(VEC_X(vec, idx));
  VEC_X(vec, idx) = DBL_Create(val);

  return STATUS_SUCCESS;
}

/*! FUNCTION:  VECTOR_DBL_Insert()
 *  SYNOPSIS:  Overwrite data from <vec> at the <idx> position in array to <val>. Deletes present value.
 */
STATUS_FLAG
VECTOR_DBL_Insert(VECTOR_DBL* vec,
                  int idx,
                  DBL val) {
/* if debugging, do edgebound checks */
#if SAFE
  if (idx >= vec->N || idx < 0) {
    fprintf(stderr, "ERROR: VECTOR_DBL access out-of-bounds.\n");
    fprintf(stderr, "dim: (%ld/%ld), access: %d\n", vec->N, vec->Nalloc, idx);
    ERRORCHECK_exit(EXIT_FAILURE);
  }
#endif

  VEC_X(vec, idx) = DBL_Destroy(VEC_X(vec, idx));
  VEC_X(vec, idx) = DBL_Create(val);

  return STATUS_SUCCESS;
}

/*! FUNCTION:  VECTOR_DBL_Delete()
 *  SYNOPSIS:  Overwrite data from <vec> at the <idx> position in array to <val>. Deletes present value.
 */
STATUS_FLAG
VECTOR_DBL_Delete(VECTOR_DBL* vec,
                  int idx) {
/* if debugging, do edgebound checks */
#if SAFE
  if (idx >= vec->N || idx < 0) {
    fprintf(stderr, "ERROR: VECTOR_DBL access out-of-bounds.\n");
    fprintf(stderr, "dim: (%ld/%ld), access: %d\n", vec->N, vec->Nalloc, idx);
    ERRORCHECK_exit(EXIT_FAILURE);
  }
#endif

  int N = VECTOR_DBL_GetSize(vec);
  VEC_X(vec, idx) = DBL_Destroy(VEC_X(vec, idx));
  VEC_X(vec, idx) = VEC_X(vec, N - 1);
  // VEC_X( vec, N - 1 )   = DBL_Empty();
  vec->N -= 1;

  return STATUS_SUCCESS;
}

/*! FUNCTION:  VECTOR_DBL_Push()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array.
 *             WARNING: Does not handle resizing or check for out-of-bounds.
 *                      For that, use Pushback().
 */
inline STATUS_FLAG
VECTOR_DBL_Push(VECTOR_DBL* vec,
                DBL val) {
  /* NOTE: This push() creates another copy of the data to store in vector (in the case of dynamically allocated data) */
  VEC_X(vec, vec->N) = DBL_Create(val);
  vec->N++;
}

/*! FUNCTION:  VECTOR_DBL_Pushback()
 *  SYNOPSIS:  Push <val> onto the end of <vec> data array,
 *             Resize array if array is full.
 */
inline STATUS_FLAG
VECTOR_DBL_Pushback(VECTOR_DBL* vec,
                    DBL val) {
  VECTOR_DBL_Push(vec, val);

  /* if array is full, resize */
  if (vec->N >= vec->Nalloc - 1) {
    VECTOR_DBL_Resize(vec, vec->N * 2);
  }
}

/*! FUNCTION:  VECTOR_DBL_Pop()
 *  SYNOPSIS:  Pop data from the end of <vec> data array, remove data, and return data.
 */
inline DBL
VECTOR_DBL_Pop(VECTOR_DBL* vec) {
  DBL data = VECTOR_DBL_Get(vec, vec->N - 1);
  vec->N -= 1;

  return data;
}

/*! FUNCTION:  VECTOR_DBL_Popback()
 *  SYNOPSIS:  Pop data from the end of <vec> data array and return data.
 *             Resize if array is less than half full.
 */
inline DBL
VECTOR_DBL_Popback(VECTOR_DBL* vec) {
  DBL data = VECTOR_DBL_Pop(vec);

  /* if array is less than half used, resize */
  if (vec->N < vec->Nalloc / 2) {
    VECTOR_DBL_Resize(vec, vec->N / 2);
  }

  return data;
}

/*! FUNCTION:  VECTOR_DBL_Append()
 *  SYNOPSIS:  Push <append> data array of length <L> onto the end of <vec> data array.
 */
inline STATUS_FLAG
VECTOR_DBL_Append(VECTOR_DBL* vec,
                  DBL* append,
                  size_t L) {
  size_t N_new;
  N_new = vec->N + L;

  /* resize array */
  if (vec->Nalloc < N_new) {
    VECTOR_DBL_Resize(vec, N_new);
  }

  /* copy data over */
  for (int i = 0; i < L; i++) {
    VECTOR_DBL_Push(vec, append[i]);
  }
}

/*! FUNCTION:  VECTOR_DBL_GetSize()
 *  SYNOPSIS:  Get utilized length of <vec>.
 */
inline int
VECTOR_DBL_GetSize(VECTOR_DBL* vec) {
  return vec->N;
}

/*! FUNCTION:  VECTOR_DBL_SetSize()
 *  SYNOPSIS:  Set utilized length of <vec>.
 *             Will allocate memory if necessary.
 */
inline STATUS_FLAG
VECTOR_DBL_SetSize(VECTOR_DBL* vec,
                   size_t size) {
  VECTOR_DBL_GrowTo(vec, size);
  vec->N = size;
}

/*! FUNCTION:  VECTOR_DBL_GetSizeAlloc()
 *  SYNOPSIS:  Get allocated length of <vec>.
 */
inline int
VECTOR_DBL_GetSizeAlloc(VECTOR_DBL* vec) {
  return vec->Nalloc;
}

/*! FUNCTION:  VECTOR_DBL_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array.
 *             Returns index of the first occurance found.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of occurance of <val>.
 *             Return -1 if <val> is not found.
 */
int VECTOR_DBL_Search(VECTOR_DBL* vec,
                      DBL val) {
  int N = VECTOR_DBL_GetSize(vec);
  int idx = (N / 2);
  int found = -1;

  for (int i = N / 4; i >= 1; i /= 2) {
    int cmp = DBL_Compare(val, vec->data[i]);

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

/*! FUNCTION:  VECTOR_DBL_Search_First()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array.
 *             Returns index of the first occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of first occurance of <val>.
 *             Return -1 if <val> is not found.
 */
int VECTOR_DBL_Search_First(VECTOR_DBL* vec,
                            DBL val) {
  int N = VECTOR_DBL_GetSize(vec);
  int idx = (N / 2);
  int found = -1;

  for (int i = N / 4; i >= 1; i /= 2) {
    int cmp = DBL_Compare(val, vec->data[i]);

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

/*! FUNCTION:  VECTOR_DBL_Search_Last()
 *  SYNOPSIS:  Binary search of <vec> to find <val> in data array.
 *             Returns index of the last occurance.
 *             Assumes <vec> has been sorted ascending.
 *  RETURN:    Returns index of last occurance of <val>.
 *             Return -1 if <val> is not found.
 */
int VECTOR_DBL_Search_Last(VECTOR_DBL* vec,
                           DBL val) {
  int N = VECTOR_DBL_GetSize(vec);
  int idx = (N / 2);
  int found = -1;

  for (int i = N / 4; i >= 1; i /= 2) {
    int cmp = DBL_Compare(val, vec->data[i]);

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

/*! FUNCTION:  VECTOR_DBL_Compare()
 *  SYNOPSIS:  Compare <vec_A> and <vec_B>.
 *  RETURN:    0 for equality,
 *             pos if <vec_A> > <vec_B>,
 *             neg if <vec_A> < <vec_B>.
 */
int VECTOR_DBL_Compare(VECTOR_DBL* vec_A,
                       VECTOR_DBL* vec_B) {
  for (int i = 0; i < vec_A->N; i++) {
    if (DBL_Compare(vec_A->data[i], vec_B->data[i]) != 0) {
      return DBL_Compare(vec_A->data[i], vec_B->data[i]);
    }
  }
  return 0;
}

/*! FUNCTION:  VECTOR_DBL_Sort()
 *  SYNOPSIS:  Sorts <vec> data array in ascending order. In-place.
 */
STATUS_FLAG
VECTOR_DBL_Sort(VECTOR_DBL* vec) {
  int N = VECTOR_DBL_GetSize(vec);
  qsort(vec->data, N, sizeof(DBL), DBL_CompareTo);

  // VECTOR_DBL_Sort_Sub( vec, 0, N );

#if DEBUG
  {
    for (int i = 0; i < N - 1; i++) {
      DBL cur = vec->data[i];
      DBL nxt = vec->data[i + 1];
      char s_cur[50];
      char s_nxt[50];
      int cmp = DBL_Compare(cur, nxt);
      if ((cmp <= 0) == false) {
        fprintf(stderr, "ERROR: bad sort. %d, %d v %d: %s vs %s\n",
                cmp, i, i + 1, DBL_ToString(cur, s_cur), DBL_ToString(nxt, s_nxt));
        ERRORCHECK_exit(EXIT_FAILURE);
      }
    }
  }
#endif
}

/*! FUNCTION:  VECTOR_DBL_Sort_Sub()
 *  SYNOPSIS:  Sorts subarray of <vec> data in ascending order on range (beg,end].
 *             Uses quicksort until length of subarray falls below threshold, then selection sort.
 */
STATUS_FLAG
VECTOR_DBL_Sort_Sub(VECTOR_DBL* vec,
                    int beg,
                    int end) {
  const int begin_select_sort = INT_MAX;
  int N = VECTOR_DBL_GetSize(vec);
  if (N <= 1) {
    return STATUS_SUCCESS;
  }

  int size = end - beg;
  /* run selection sort if below threshold */
  if (size <= begin_select_sort) {
    VECTOR_DBL_Sort_Sub_Selectsort(vec, beg, end);
  }
  /* otherwise run quiksort */
  else {
    VECTOR_DBL_Sort_Sub_Quicksort(vec, beg, end);
  }

  return STATUS_SUCCESS;
}

/*! FUNCTION:  VECTOR_DBL_Sort_Sub_Selectsort()
 *  SYNOPSIS:  Selection Sorts subarray of <vec> data in ascending order on range (beg,end].
 */
STATUS_FLAG
VECTOR_DBL_Sort_Sub_Selectsort(VECTOR_DBL* vec,
                               int beg,
                               int end) {
  /* find the minimum element of remaining unsorted list */
  for (int i = beg; i < end; i++) {
    /* initial minimum value found */
    int min_idx = i;
    DBL min_val = vec->data[i];
    for (int j = i + 1; j < end; j++) {
      /* if new minimum found, update value and index */
      int cmp = DBL_Compare(min_val, vec->data[j]);
      if (cmp > 0) {
        min_idx = j;
        min_val = vec->data[j];
      }
    }
    /* swap new minimum to left-most position */
    VECTOR_DBL_Swap(vec, i, min_idx);
  }
}

/*! FUNCTION:  VECTOR_DBL_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Quick Sorts subarray of <vec> data in ascending order on range (beg,end].
 */
STATUS_FLAG
VECTOR_DBL_Sort_Sub_Quicksort(VECTOR_DBL* vec,
                              int beg,
                              int end) {
  /* partition pointers */
  int r_idx = beg + 1;
  int l_idx = end - 1;
  DBL* rhs = &(vec->data[beg + 1]);
  DBL* lhs = &(vec->data[end - 1]);

  /* select random pivot value */
  int range = end - beg;
  int pivot_idx = RNG_INT_Range(beg, end);
  DBL pivot_val = vec->data[pivot_idx];
  VECTOR_DBL_Swap(vec, pivot_idx, beg);

  /* partition on pivot */
  while (l_idx <= r_idx) {
    /* find next right partition element that is less than pivot element */
    while ((l_idx <= r_idx) && (DBL_Compare(pivot_val, vec->data[r_idx]) < 0)) {
      r_idx--;
    }
    /* find next left partition element that is greater than pivot element */
    while ((l_idx <= r_idx) && (DBL_Compare(pivot_val, vec->data[l_idx]) >= 0)) {
      l_idx++;
    }
    /* if left and right index have not crossed, then swap elements */
    if (l_idx <= r_idx) {
      VECTOR_DBL_Swap(vec, l_idx, r_idx);
    }
  }
  /* move partition element to barrier between left and right index */
  VECTOR_DBL_Swap(vec, beg, r_idx);
  /* sort both partitions (omit partition element) */
  VECTOR_DBL_Sort_Sub(vec, beg, r_idx);
  VECTOR_DBL_Sort_Sub(vec, r_idx + 1, end);
}

/*! FUNCTION:  VECTOR_DBL_Op()
 *  SYNOPSIS:  Perform element-wise unary operation <op> to each cell in <vec_in> and puts it in <vec_out>.
 *             Returns a pointer to <vec_out>.
 *             <vec_out> will be resized to size of <vec_in>.
 *             <vec_in> and <vec_out> can be the same vector. If <vec_out> is NULL, new vector will be created.
 */
inline VECTOR_DBL*
VECTOR_DBL_Op(VECTOR_DBL* vec_out, /* input vector */
              VECTOR_DBL* vec_in,  /* output vector (can be input vector) */
              DBL (*op)(DBL data)) /* unary operation */
{
  int N = VECTOR_DBL_GetSize(vec_in);
  if (vec_out == NULL) {
    VECTOR_DBL_Create_by_Size(N);
  }
  VECTOR_DBL_SetSize(vec_out, N);

  for (int i = 0; i < N; i++) {
    *VECTOR_DBL_GetX(vec_out, i) = op(*VECTOR_DBL_GetX(vec_in, i));
  }

  return vec_out;
}

/*! FUNCTION:  VECTOR_DBL_Op()
 *  SYNOPSIS:  Perform element-wise binary operation <op>(DBL data_1, DBL data_2) to each cell in <vec_in_1, vec_in_2> and puts it in <vec_out>.
 *             <vec_in_1> and <vec_in_2> must be the same size.
 *             Returns a pointer to <vec_out>.
 *             <vec_out> will be resized to size of <vec_in>.
 *             <vec_in> and <vec_out> can be the same vector. If <vec_out> is NULL, new vector will be created.
 */
inline VECTOR_DBL*
VECTOR_DBL_BinOp(VECTOR_DBL* vec_out,               /* output vector */
                 VECTOR_DBL* vec_in_1,              /* first input vector */
                 VECTOR_DBL* vec_in_2,              /* second input vector */
                 DBL (*op)(DBL data_1, DBL data_2)) /* binary operation */
{
  int N = VECTOR_DBL_GetSize(vec_in_1);
  int N_2 = VECTOR_DBL_GetSize(vec_in_2);
  if (N != N_2) {
    fprintf(stderr, "ERROR: <vec_in_1> and <vec_in_2> are not the same dimensions.\n");
    ERRORCHECK_exit(EXIT_FAILURE);
  }
  if (vec_out == NULL) {
    VECTOR_DBL_Create_by_Size(N);
  }
  VECTOR_DBL_SetSize(vec_out, N);

  for (int i = 0; i < N; i++) {
    *VECTOR_DBL_GetX(vec_out, i) = op(*VECTOR_DBL_GetX(vec_in_1, i), *VECTOR_DBL_GetX(vec_in_2, i));
  }

  return vec_out;
}

/*! FUNCTION:  VECTOR_DBL_Swap()
 *  SYNOPSIS:  Swaps the values of <vec> at indexes <i> and <j>.
 *             Warning: Only checks for Out-of-Bounds when in DEBUG.
 */
inline STATUS_FLAG
VECTOR_DBL_Swap(VECTOR_DBL* vec,
                int i,
                int j) {
  DBL swap = VECTOR_DBL_Get(vec, i);
  vec->data[i] = vec->data[j];
  vec->data[j] = swap;
}

/*! FUNCTION:  VECTOR_DBL_Reverse()
 *  SYNOPSIS:  Reverse the ordering of array.
 */
inline STATUS_FLAG VECTOR_DBL_Reverse(VECTOR_DBL* vec) {
  int N = VECTOR_DBL_GetSize(vec);

  for (int i = 0; i < (N / 2); i++) {
    VECTOR_DBL_Swap(vec, i, (N - 1) - i);
  }
}

/*! FUNCTION:  VECTOR_DBL_Dump()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer. Non-optimized.
 */
STATUS_FLAG
VECTOR_DBL_Dump(VECTOR_DBL* vec,
                FILE* fp) {
  VECTOR_DBL_Dump_byOpt(vec, "\n", "VECTOR BOUND", fp);
}

/*! FUNCTION:  VECTOR_DBL_Dump_byOpt()
 *  SYNOPSIS:  Output <vec> to <fp> file pointer.
 */
STATUS_FLAG
VECTOR_DBL_Dump_byOpt(VECTOR_DBL* vec,
                      STR delim,
                      STR header,
                      FILE* fp) {
  /* stringification of template object */
  char s[50];
  const char* pad = " ";

  fprintf(fp, "%s: ", header);
  fprintf(fp, "[ ");
  for (int i = 0; i < vec->N; i++) {
    fprintf(fp, "%s%s%s", DBL_ToString(vec->data[i], s), delim, pad);
  }
  if (vec->N >= 1) {
    fprintf(fp, "%s%s", DBL_ToString(vec->data[vec->N - 1], s), pad);
  }
  fprintf(fp, "]\n");
}

/*! FUNCTION:  VECTOR_DBL_UnitTest()
 *  SYNOPSIS:  Perform unit test for VECTOR_DBL.
 */
STATUS_FLAG
VECTOR_DBL_UnitTest() {
  VECTOR_DBL* vec = VECTOR_DBL_Create();

  VECTOR_DBL_Destroy(vec);
}
