/*******************************************************************************
 *  @file edgebound.c
 *  @brief Edgebounds Datatype
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/* local imports */
#include "../structs.h"
#include "edgebound.h"

/*
 *  FUNCTION: edgebounds_Create()
 *  SYNOPSIS: Create new edgebounds object and return pointer.
 *
 *  PURPOSE:   
 *
 *  ARGS:      
 *
 *  RETURN:    <edg>      Edgebounds Object
 */
EDGEBOUNDS *edgebounds_Create()
{
   EDGEBOUNDS *edg;
   const int min_size = 16;
   edg = (EDGEBOUNDS *)malloc(sizeof(EDGEBOUNDS));
   edg->N = 0;
   edg->size = min_size;
   edg->bounds = (BOUND *)malloc(min_size * sizeof(BOUND));
   return edg;
}

/*
 *  FUNCTION: edgebounds_Init()
 *  SYNOPSIS: Initialize new edgebounds object pointer.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>      Edgebounds Object
 *
 *  RETURN:
 */
void edgebounds_Init(EDGEBOUNDS **edg)
{
   const int min_size = 16;
   (*edg) = (EDGEBOUNDS *)malloc(sizeof(EDGEBOUNDS));
   (*edg)->N = 0;
   (*edg)->size = min_size;
   (*edg)->bounds = (BOUND *)malloc(min_size * sizeof(BOUND));
}

/*
 *  FUNCTION: edgebounds_Destroy()
 *  SYNOPSIS: Destroy edgebounds object.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>      Edgebounds Object
 *
 *  RETURN:
 */
void edgebounds_Destroy(EDGEBOUNDS *edg)
{
   free(edg->bounds);
   free(edg);
}


/*
 *  FUNCTION: edgebounds_Push()
 *  SYNOPSIS: Add bound to Edgebound list.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>       Edgebounds,
 *             <bnd>       Bound
 *
 *  RETURN:
 */
void edgebounds_Add(EDGEBOUNDS *edg,
                    BOUND bnd)
{
   edg->bounds[edg->N].diag = bnd.diag;
   edg->bounds[edg->N].lb = bnd.lb;
   edg->bounds[edg->N].rb = bnd.rb;

   /* resize if necessary */
   edg->N += 1;
   if (edg->N >= edg->size) 
      edgebounds_Resize(edg);
}


/*
 *  FUNCTION: edgebounds_Insert()
 *  SYNOPSIS: Insert bound into i-index of Edgebound list.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>       Edgebounds,
 *             <bnd>       Bound,
 *             <i>         int Index
 *  RETURN:
 */
void edgebounds_Insert(EDGEBOUNDS *edg,
                       BOUND bnd,
                       int i)
{
   edg->bounds[i].diag = bnd.diag;
   edg->bounds[i].lb = bnd.lb;
   edg->bounds[i].rb = bnd.rb;
}


/*
 *  FUNCTION: edgebounds_Delete()
 *  SYNOPSIS: Delete bound at i-index and fill from end of list.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>       Edgebounds,
 *             <bnd>       Bound,
 *             <i>         int Index
 *  RETURN:
 */
void edgebounds_Delete(EDGEBOUNDS *edg,
                       BOUND bnd,
                       int i)
{
   int N = edg->N;
   edg->bounds[i].diag = edg->bounds[N].diag;
   edg->bounds[i].lb = edg->bounds[N].lb;
   edg->bounds[i].rb = edg->bounds[N].rb;
   edg->N -= 1;
}

/*
 *  FUNCTION: edgebounds_Resize()
 *  SYNOPSIS: Resize number of bounds in edgebound object.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>      Edgebounds Object
 *
 *  RETURN:
 */
void edgebounds_Resize(EDGEBOUNDS *edg)
{
   const int growth_factor = 2;
   edg->size *= growth_factor;
   edg->bounds = (BOUND *)realloc(edg->bounds, edg->size * sizeof(BOUND));
}


/*
 *  FUNCTION: edgebounds_Reverse()
 *  SYNOPSIS: Reverse order of edgebound list.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>      Edgebounds Object
 *
 *  RETURN:
 */
void edgebounds_Reverse(EDGEBOUNDS *edg)
{
   BOUND tmp;
   for (int i = 0; i <= (edg->N / 2); ++i)
   {
      tmp.diag = edg->bounds[i].diag;
      tmp.lb = edg->bounds[i].lb;
      tmp.rb = edg->bounds[i].rb;

      edg->bounds[i].diag = edg->bounds[edg->N-i].diag;
      edg->bounds[i].lb = edg->bounds[edg->N-i].lb;
      edg->bounds[i].rb = edg->bounds[edg->N-i].rb;

      edg->bounds[edg->N-i].diag = tmp.diag;
      edg->bounds[edg->N-i].lb = tmp.lb;
      edg->bounds[edg->N-i].rb = tmp.rb;
   }
}


/*
 *  FUNCTION: edgebounds_Print()
 *  SYNOPSIS: Print EDGEBOUND object.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>      Edgebounds Object
 *
 *  RETURN:
 */
void edgebounds_Print(EDGEBOUNDS *edg)
{
   printf("printing edgebounds...\n");
   printf("N: %d, Nalloc: %d\n", edg->N, edg->size);
   for (unsigned int i = 0; i < edg->N; ++i)
   {
      printf("[%d] ", i);
      bound_Print(edg->bounds[i]);
   }
}


/*
 *  FUNCTION: bound_Print()
 *  SYNOPSIS: Print BOUND object.
 *
 *  PURPOSE:
 *
 *  ARGS:      <bnd>      Bounds Object
 *
 *  RETURN:
 */
void bound_Print(BOUND bnd)
{
   printf("{ d: %d, lb: %d, rb: %d }\n", bnd.diag, bnd.lb, bnd.rb);
}


/*
 *  FUNCTION: edgebounds_Save()
 *  SYNOPSIS: Save edgebound printout to file.
 *
 *  PURPOSE:
 *
 *  ARGS:      <bnd>      Bounds Object
 *             <f>        Filename
 *
 *  RETURN:
 */
void edgebounds_Save(EDGEBOUNDS *edg,
                      const char *_filename_)
{
   FILE *fp;
   fp = fopen(_filename_, "w");

   for (unsigned int i = 0; i < edg->N; ++i)
   {
      fprintf(fp, "[%d] x: %d\t y: (%d, %d)\n", i, edg->bounds[i].diag, edg->bounds[i].lb, edg->bounds[i].rb);
   }
   fclose(fp);
}


/*
 *  FUNCTION: bounds_Compare()
 *  SYNOPSIS: Compare two Bounds, first by diag, then by lb, then by rb
 *
 *  PURPOSE:
 *
 *  ARGS:      <a>        Bound,
 *             <b>        Bound
 *
 *  RETURN:    1 if (a > b), 0 if equal, -1 if (a < b)
 */
int bounds_Compare(BOUND a, BOUND b)
{
   if (a.diag > b.diag) {
      return 1;
   } else 
   if (a.diag < b.diag) {
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