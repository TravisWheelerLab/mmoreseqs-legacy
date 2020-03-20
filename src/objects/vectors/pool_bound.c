/*******************************************************************************
 *  FILE:      POOL_BOUND.c
 *  SYNOPSIS:  OBJECT POOL of BOUND data type.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

/* local imports */
#include "../structs.h"
#include "vectors/vector_char.h"
#include "vectors/vector_int.h"
#include "vectors/vector_bound.h"

/* header */
#include "pool_bound.h"

/* constructor */
POOL_BOUND* POOL_BOUND_Create()
{
   const int min_size = 8;
   POOL_BOUND* pool;

   pool = (POOL_BOUND*) malloc( sizeof(POOL_BOUND) );
   if (pool == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc POOL_BOUND object.\n");
      exit(EXIT_FAILURE);
   }

   pool->N       = 0;
   pool->Nalloc  = min_size;
   pool->isOcc   = VECTOR_CHAR_Create();
   pool->heads   = VECTOR_INT_Create();
   pool->bounds  = VECTOR_BOUND_Create();

   return pool;
}


/* destructor */
void POOL_BOUND_Destroy(POOL_BOUND* pool)
{
   VECTOR_CHAR_Destroy(pool->isOcc);
   VECTOR_INT_Destroy(pool->heads);
   VECTOR_BOUND_Destroy(pool->bounds);
   free(pool);
}