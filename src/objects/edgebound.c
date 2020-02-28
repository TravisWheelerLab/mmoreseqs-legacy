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
#include "bound.h"

/* header */
#include "edgebound.h"

/*
 *  FUNCTION: EDGEBOUNDS_Create()
 *  SYNOPSIS: Create new EDGEBOUNDS object and returns pointer.
 *
 *  ARGS:      None.
 *
 *  RETURN:    <edg>      Edgebounds Object
 */
EDGEBOUNDS* EDGEBOUNDS_Create()
{
   const int min_size = 8;
   EDGEBOUNDS *edg;
   edg = (EDGEBOUNDS*)malloc(sizeof(EDGEBOUNDS));
   if (edg == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc for EDGEBOUNDS.\n");
   }

   edg->N      = 0;
   edg->Nalloc = min_size;
   edg->ids    = VECTOR_INT_Create();
   edg->heads  = VECTOR_INT_Create();

   edg->bounds = (BOUND*)malloc(min_size * sizeof(BOUND));
   if (edg == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc BOUNDS for EDGEBOUNDS.\n");
   }
   return edg;
}


/*
 *  FUNCTION: EDGEBOUNDS_Create()
 *  SYNOPSIS: Create new EDGEBOUNDS object and returns pointer.
 *
 *  ARGS:      None.
 *
 *  RETURN:    <edg>      Edgebounds Object
 */
void EDGEBOUNDS_Init(EDGEBOUNDS **edg)
{
   *edg = EDGEBOUNDS_Create();
}

/*
 *  FUNCTION: EDGEBOUNDS_Destroy()
 *  SYNOPSIS: Frees all memory from EDGEBOUNDS object.
 *
 *  ARGS:     <edg>      Edgebounds Object
 *
 *  RETURN:   None.
 */
void EDGEBOUNDS_Destroy(EDGEBOUNDS*  edg)
{
   free(edg->bounds);
   free(edg);
}


/*
 *  FUNCTION: EDGEBOUNDS_Pushback()
 *  SYNOPSIS: Add bound to Edgebound list.
 *
 *  ARGS:      <edg>       Edgebounds,
 *             <bnd>       Bound
 *
 *  RETURN:
 */
void EDGEBOUNDS_Pushback(EDGEBOUNDS*  edg,
                         BOUND        bnd)
{
   /* resize if necessary */
   if (edg->N >= edg->Nalloc - 1) 
      EDGEBOUNDS_Resize(edg);

   edg->bounds[edg->N] = bnd;
   edg->N++;
}


/*
 *  FUNCTION: EDGEBOUNDS_Pushback_Head()
 *  SYNOPSIS: Add head index and row/diag id to Head list.
 *
 *  ARGS:      <edg>       Edgebounds,
 *             <id>        id for row/diag,
 *             <head>      head index for id in row/diag
 *
 *  RETURN:
 */
void EDGEBOUNDS_Pushback_Head(EDGEBOUNDS* edg,
                              int         id,
                              int         head)
{
   VECTOR_INT_Pushback(edg->ids, id);
   VECTOR_INT_Pushback(edg->heads, head);
}


/*
 *  FUNCTION: EDGEBOUNDS_Insert()
 *  SYNOPSIS: Insert/Overwrite bound into i-index of Edgebound list.
 *
 *  ARGS:      <edg>       Edgebounds,
 *             <bnd>       Bound,
 *             <i>         int Index
 *  RETURN:
 */
void EDGEBOUNDS_Insert(EDGEBOUNDS*  edg,
                       BOUND        bnd,
                       int          i)
{
   edg->bounds[i] = bnd;
}


/*
 *  FUNCTION: EDGEBOUNDS_Delete()
 *  SYNOPSIS: Delete bound at i-index and fill from end of list.
 *
 *  ARGS:      <edg>       Edgebounds,
 *             <bnd>       Bound,
 *             <i>         int Index
 *  RETURN:
 */
void EDGEBOUNDS_Delete(EDGEBOUNDS*  edg,
                       BOUND        bnd,
                       int          i)
{
   int N = edg->N;
   edg->bounds[i].id = edg->bounds[N].id;
   edg->bounds[i].lb = edg->bounds[N].lb;
   edg->bounds[i].rb = edg->bounds[N].rb;
   edg->N -= 1;
}

/*
 *  FUNCTION: EDGEBOUNDS_Resize()
 *  SYNOPSIS: Resize number of bounds in edgebound object.
 *
 *  ARGS:      <edg>      Edgebounds Object
 *
 *  RETURN:
 */
void EDGEBOUNDS_Resize(EDGEBOUNDS *edg)
{
   const int growth_factor = 2;
   edg->Nalloc *= growth_factor;
   edg->bounds = (BOUND *)realloc(edg->bounds, edg->Nalloc * sizeof(BOUND));
}

/*
 *  FUNCTION: EDGEBOUNDS_Reverse()
 *  SYNOPSIS: Reverse order of edgebound list.
 *
 *  ARGS:      <edg>      Edgebounds Object
 *
 *  RETURN:
 */
void EDGEBOUNDS_Reverse(EDGEBOUNDS *edg)
{
   BOUND tmp;
   for (int i = 0; i <= (edg->N / 2); ++i)
   {
      tmp.id = edg->bounds[i].id;
      tmp.lb = edg->bounds[i].lb;
      tmp.rb = edg->bounds[i].rb;

      edg->bounds[i].id = edg->bounds[edg->N-i].id;
      edg->bounds[i].lb = edg->bounds[edg->N-i].lb;
      edg->bounds[i].rb = edg->bounds[edg->N-i].rb;

      edg->bounds[edg->N-i].id = tmp.id;
      edg->bounds[edg->N-i].lb = tmp.lb;
      edg->bounds[edg->N-i].rb = tmp.rb;
   }
}

/*
 *  FUNCTION: EDGEBOUNDS_SetHeads()
 *  SYNOPSIS: Traverse Bounds and Find Heads of Each Row (Boundaries)
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>      Edgebounds Object
 *
 *  RETURN:
 */
void EDGEBOUNDS_SetHeads(EDGEBOUNDS *edg)
{
   int id;
   id = edg->bounds[0].id;
   VECTOR_INT_Pushback(edg->ids, id);
   VECTOR_INT_Pushback(edg->heads, 0);

   for (int i = 1; i < edg->N; i++) {
      if (edg->bounds[i-1].id != edg->bounds[i].id) {
         id = edg->bounds[i].id;
         VECTOR_INT_Pushback(edg->ids, id);
         VECTOR_INT_Pushback(edg->heads, i);
      }
   }
}

/*
 *  FUNCTION: EDGEBOUNDS_Print()
 *  SYNOPSIS: Print EDGEBOUND object.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>      Edgebounds Object
 *
 *  RETURN:
 */
void EDGEBOUNDS_Dump(EDGEBOUNDS* edg,
                     FILE*       fp)
{
   /* test for bad file pointer */
   if (fp == NULL) {
      const char* obj_name = "EDGEBOUNDS";
      fprintf(stderr, "ERROR: Bad FILE POINTER for printing %s.\n", obj_name);
      exit(EXIT_FAILURE);
      return;
   }
   printf("printing edgebound...\n");
   fprintf(fp, "N: %d, Nalloc: %d\n", edg->N, edg->Nalloc);
   for (unsigned int i = 0; i < edg->N; ++i)
   {
      BOUND bnd = edg->bounds[i];
      fprintf(fp, "[%d] ", i);
      fprintf(fp, "{ id: %d, lb: %d, rb: %d }\n", bnd.id, bnd.lb, bnd.rb);
   }
}

/*
 *  FUNCTION: EDGEBOUNDS_Dump()
 *  SYNOPSIS: Save edgebound printout to file.
 *
 *  ARGS:      <bnd>      Bounds Object
 *             <f>        Filename
 *
 *  RETURN:
 */
void EDGEBOUNDS_Save(EDGEBOUNDS*  edg,
                     const char*  _filename_)
{
   FILE *fp;
   fp = fopen(_filename_, "w");
   EDGEBOUNDS_Dump(edg, fp);
   fclose(fp);
}



/*
 *  FUNCTION: EDGEBOUNDS_Count()
 *  SYNOPSIS: Count the number of cells in edgebound.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>      Edgebounds Object
 *
 *  RETURN:
 */
int EDGEBOUNDS_Count(EDGEBOUNDS *edg)
{
   int sum = 0;
   for (int i = 0; i < edg->N; i++)
   {
      sum += edg->bounds[i].rb - edg->bounds[i].lb;
   }
   return sum;
}