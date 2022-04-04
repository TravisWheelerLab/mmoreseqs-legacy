/*******************************************************************************
 *  - FILE:   matrix_2d.h
 *  - DESC:    MATRIX_2D Float object.
 *******************************************************************************/

#ifndef _MATRIX_2D_H
#define _MATRIX_2D_H

/*! FUNCTION:  MATRIX_2D_Create()
 *  SYNOPSIS:  Constructor: Allocates memory for MATRIX_2D.
 *
 *  RETURN:    Return pointer to new MATRIX_2D object.
 */
MATRIX_2D* MATRIX_2D_Create(int R, int C);

/*! FUNCTION:  MATRIX_2D_Create_Clean()
 *  SYNOPSIS:  Constructor: Allocates memory for MATRIX_2D. Sets all data to
 * -INF.
 *
 *  RETURN:    Return pointer to new MATRIX_2D object.
 */
MATRIX_2D* MATRIX_2D_Create_Clean(int R, int C);

/*! FUNCTION:  MATRIX_2D_Destroy()
 *  SYNOPSIS:  Destructor: Frees memory for MATRIX_2D.
 *
 *  RETURN:    Return pointer to new MATRIX_2D object.
 */
MATRIX_2D* MATRIX_2D_Destroy(MATRIX_2D* mx);

/*! FUNCTION:  MATRIX_2D_Copy()
 *  SYNOPSIS:  Deep Copy: Makes deep copy of <src> to <dest>.  Creates new
 * MATRIX_2D object if <dest> is NULL.
 *
 *  RETURN:    Return pointer to <dest>.
 */
MATRIX_2D* MATRIX_2D_Copy(MATRIX_2D* dest, const MATRIX_2D* src);

/*  FUNCTION:  MATRIX_2D_Add()
 *  SYNOPSIS:  Takes sum of <mx_A> + <mx_B>.  Result stored in <mx_res>.
 */
int MATRIX_2D_Add(MATRIX_2D* mx_A, MATRIX_2D* mx_B, MATRIX_2D* mx_res);

/*! FUNCTION:
 *  SYNOPSIS:
 *
 *  RETURN:
 */
int MATRIX_2D_Fill(MATRIX_2D* mx, float val);

/* fill MATRIX_2D with -INF */
/*! FUNCTION:
 *  SYNOPSIS:
 *
 *  RETURN:
 */
int MATRIX_2D_Clean(MATRIX_2D* mx);

/* check that all cells are filled with given value */
/*! FUNCTION:
 *  SYNOPSIS:
 *
 *  RETURN:
 */
int MATRIX_2D_Check_Value(MATRIX_2D* mx, float val);

/* Check MATRIX_2D data filled with -INF */
/*! FUNCTION:
 *  SYNOPSIS:
 *
 *  RETURN:
 */
int MATRIX_2D_Check_Clean(MATRIX_2D* mx);

/* getter for index */
/*! FUNCTION:
 *  SYNOPSIS:
 *
 *  RETURN:
 */
float* MATRIX_2D_Get(MATRIX_2D* mx, int i, int j);

/* getter pointer for index in MATRIX (input in 1D-coords) */
/*! FUNCTION:
 *  SYNOPSIS:
 *
 *  RETURN:
 */
float* MATRIX_2D_Get1D(MATRIX_2D* mx, int n);

/* convert 2D-coords to 1D-coords */
/*! FUNCTION:
 *  SYNOPSIS:
 *
 *  RETURN:
 */
int MATRIX_2D_to_1D(MATRIX_2D* mx, int i, int j);

/* reuse matrix by resizing only if new matrix requires more memory */
/*! FUNCTION:
 *  SYNOPSIS:
 *
 *  RETURN:
 */
int MATRIX_2D_Reuse(MATRIX_2D* mx, int R, int C);

/* reuse MATRIX_2D by resizing only if new matrix requires more memory.  All new
 * matrix values are cleaned. */
/*! FUNCTION:
 *  SYNOPSIS:
 *
 *  RETURN:
 */
int MATRIX_2D_Reuse_Clean(MATRIX_2D* mx, int R, int C);

/* resize matrix  */
/*! FUNCTION:
 *  SYNOPSIS:
 *
 *  RETURN:
 */
int MATRIX_2D_Resize(MATRIX_2D* mx, int R, int C);

/* Outputs MATRIX_2D out to FILE POINTER */
/*! FUNCTION:
 *  SYNOPSIS:
 *
 *  RETURN:
 */
int MATRIX_2D_Dump(MATRIX_2D* mx, FILE* fp);

/* Save MATRIX_2D to FILE by FILENAME */
/*! FUNCTION:
 *  SYNOPSIS:
 *
 *  RETURN:
 */
int MATRIX_2D_Save(MATRIX_2D* mx, char* filename);

/* Compare two MATRIX_2D */
/*! FUNCTION:
 *  SYNOPSIS:
 *
 *  RETURN:
 */
int MATRIX_2D_Compare(MATRIX_2D* mx_A, MATRIX_2D* mx_B);

/*
 *  FUNCTION:  MATRIX_2D_Diff()
 *  SYNOPSIS:  Takes difference of <mx_A> - <mx_B>.  Result stored in <mx_diff>.
 */
/*! FUNCTION:
 *  SYNOPSIS:
 *
 *  RETURN:
 */
int MATRIX_2D_Diff(MATRIX_2D* mx_A, MATRIX_2D* mx_B, MATRIX_2D* mx_diff);

/*! FUNCTION:  MATRIX_2D_Log()
 *  SYNOPSIS:  Performs logrithmic function log() on each cell in matrix.
 *
 *  RETURN:    None.
 */
int MATRIX_2D_Log(MATRIX_2D* mx);

/*! FUNCTION:  MATRIX_2D_Logf()
 *  SYNOPSIS:  Performs logrithmic function log() on each cell in matrix.
 *
 *  RETURN:    None.
 */
int MATRIX_2D_Logf(MATRIX_2D* mx);

/*! FUNCTION:  MATRIX_2D_Exp()
 *  SYNOPSIS:  Performs exponential function exp() on each cell in matrix.
 *
 *  RETURN:    None.
 */
int MATRIX_2D_Exp(MATRIX_2D* mx);

/* Unit Test */
/*! FUNCTION:
 *  SYNOPSIS:
 *
 *  RETURN:
 */
int MATRIX_2D_UnitTest();

#endif /* _MATRIX_2D_H */
