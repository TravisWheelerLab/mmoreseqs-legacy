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
#include "vector_bound.h"

/* header */
#include "bound.h"

/*
 *  FUNCTION: BOUND_Dump()
 *  SYNOPSIS: Output BOUND data to File Pointer.
 *
 *  ARGS:      <bnd>      Bound,
 *             <fp>       File Pointer
 *
 *  RETURN:    No Return.
 */
void BOUND_Dump(const BOUND  bnd,
                FILE*        fp)
{
   fprintf(fp, "{ id: %d, lb: %d, rb: %d }\n", bnd.id, bnd.lb, bnd.rb);
}

/*
 *  FUNCTION: BOUND_Compare()
 *  SYNOPSIS: Compare two Bounds, first by diagonal, then by left-bound, then by right-bound.
 *
 *  ARGS:      <a>        Bound,
 *             <b>        Bound
 *
 *  RETURN:    1 if (a > b), 0 if equal, -1 if (a < b)
 */
int BOUND_Compare(BOUND a, 
                  BOUND b)
{
   if (a.id > b.id) {
      return 1;
   } else 
   if (a.id < b.id) {
      return -1;
   }

   if (a.lb > b.lb) {
      return 1;
   } else
   if (a.lb < b.lb) {
      return -1;
   }

   if (a.rb > b.rb) {
      return 1;
   } else
   if (a.rb < b.rb) {
      return -1;
   }
   
   return 0;
}