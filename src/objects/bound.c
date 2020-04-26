/*******************************************************************************
 *  FILE:      edgebound.c
 *  PURPOSE:   EDGEBOUNDS Object
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/* local imports */
#include "structs.h"
#include "utilities.h"
#include "objects.h"

/* header */
#include "bound.h"

/*
 *  FUNCTION:  BOUND_Dump()
 *  SYNOPSIS:  Output BOUND data to File Pointer (checks for valid pointer).
 */
void BOUND_Dump(const BOUND*     bnd,
                FILE*            fp)
{
   fprintf(fp, "{ id: %d, lb: %d, rb: %d }\n", bnd->id, bnd->lb, bnd->rb);
}

/*
 *  FUNCTION:  BOUND_Compare()
 *  SYNOPSIS:  Compare two Bounds.
 *   PROCESS:  [1]   Compare the id (which row or antidiag).  If equal, continue.
 *             [2]   Compare the left bound.  If equal, continue.
 *             [3]   Compare the right bound.
 *    RETURN:  1 if (a > b), 0 if equal, -1 if (a < b)
 */
int BOUND_Compare(const BOUND*   a, 
                  const BOUND*   b )
{
   if (a->id > b->id) {
      return 1;
   } else 
   if (a->id < b->id) {
      return -1;
   }

   if (a->lb > b->lb) {
      return 1;
   } else
   if (a->lb < b->lb) {
      return -1;
   }

   if (a->rb > b->rb) {
      return 1;
   } else
   if (a->rb < b->rb) {
      return -1;
   }
   
   return 0;
}