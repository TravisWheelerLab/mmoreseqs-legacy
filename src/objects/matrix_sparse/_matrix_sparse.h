/*******************************************************************************
 *  - FILE:   _matrix_sparse.h
 *  - DESC:    Sparse Matrix Data Types
 *******************************************************************************/

#ifndef _MATRIX_SPARSE_H
#define _MATRIX_SPARSE_H

/* edgebounds (used for defining sparse matrix space) */
#include "edgebound.h"
/* edgebound rows (used for reorienting edgebounds from antidiag-wise to
 * row-wise) */
#include "edgebound_rows.h"
/* methods for forming union of edgebounds and reorienting edgebounds */
#include "edgebound_merge_reorient.h"

/* sparse matrix (dependent on edgebounds) */
#include "matrix_3d_sparse.h"
/* sparse matrix builder */
#include "matrix_3d_sparse_build.h"

#endif /* _MATRIX_SPARSE_H */
