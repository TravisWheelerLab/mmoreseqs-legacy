/*******************************************************************************
 *  FILE:      edgebound.c
 *  PURPOSE:   EDGEBOUNDS Object
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/


#ifndef _BOUND_H
#define _BOUND_H

/*
 *  FUNCTION:  BOUND_Dump()
 *  SYNOPSIS:  Output BOUND data to File Pointer (checks for valid pointer).
 */
void BOUND_Dump(const BOUND*     bnd,
                FILE*            fp);

/*
 *  FUNCTION:  BOUND_Compare()
 *  SYNOPSIS:  Compare two Bounds.
 *   PROCESS:  [1]   Compare the id (which row or antidiag).  If equal, continue.
 *             [2]   Compare the left bound.  If equal, continue.
 *             [3]   Compare the right bound.
 *    RETURN:  1 if (a > b), 0 if equal, -1 if (a < b)
 */
int BOUND_Compare(const BOUND*   a, 
                  const BOUND*   b );

#endif /* _BOUND_H */