/*******************************************************************************
 *  FILE:      _matrix.h
 *  PURPOSE:   Matrix Data Types (2D, 3D, )
 *
 *  AUTHOR:    Dave Rich
 *   NOTES:
 *    - There is one tricky dependency here: MATRIX_3D_SPARSE is dependent on EDGEBOUNDS.
 *      This may mean that MATRIX_3D_SPARSE should be pulled out of /matrix/ folder.
 *    - Would like to add at least a INT matrix, if not a GENERIC matrix.  May require
 *      changes in naming conventions.
 *******************************************************************************/

#ifndef _MATRIX_H
#define _MATRIX_H

#include "matrix_2d.h"
#include "matrix_3d.h"

/* sparse matrix is dependent on EDGEBOUNDS */
#include "matrix_3d_sparse.h"

#endif /* _MATRIX_H */
