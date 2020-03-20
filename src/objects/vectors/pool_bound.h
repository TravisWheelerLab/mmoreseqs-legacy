/*******************************************************************************
 *  FILE:      POOL_BOUND.h
 *  SYNOPSIS:  OBJECT POOL of BOUND data type.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:
 *******************************************************************************/

#ifndef _POOL_BOUND_H
#define _POOL_BOUND_H

/* === OBJECTS === */

typedef struct {
   int            N;
   int            Nalloc;
   VECTOR_INT*    heads;      /* gives indices to the beginning of the objects in each row */
   VECTOR_CHAR*   isOcc;      /* tells whether a object is in the pool or not */
   VECTOR_BOUND*  bounds;     /* pool of BOUND objects */
} POOL_BOUND;

/* === FUNCTIONS === */

/* constructor */
POOL_BOUND* POOL_BOUND_Create();
/* destructor */
void POOL_BOUND_Destroy(POOL_BOUND* pool);

#endif /* _POOL_BOUND_H */