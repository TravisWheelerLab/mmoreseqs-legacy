/*******************************************************************************
 *  FILE:      alignment.c
 *  PURPOSE:   ALIGNMENT Object.
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
#include <math.h>

/* local imports */
#include "structs.h"
#include "utilities.h"
#include "objects.h"

/* header */
#include "alignment.h"

/* constructor */
ALIGNMENT* ALIGNMENT_Create()
{
   ALIGNMENT *aln       = NULL;
   const int min_size   = 256;

   aln = (ALIGNMENT*) malloc( sizeof(ALIGNMENT) );
   if (aln == NULL) {
      const char* obj_name = "ALIGNMENT";
      fprintf(stderr, "ERROR: Couldn't malloc %s: <%p>.\n", obj_name, aln);
      exit(EXIT_FAILURE);
   }

   aln->N      = 0;
   aln->Nalloc = min_size;
   aln->beg    = 0;
   aln->end    = 0;
   aln->traces = NULL;

   aln->Q      = 0;
   aln->T      = 0;

   ALIGNMENT_Resize(aln, min_size);

   return aln;
}

/* destructor */
void ALIGNMENT_Destroy(ALIGNMENT* aln)
{
   if (aln == NULL) return;
   
   free(aln->traces);
   free(aln);
   aln = NULL;
}

/* reuse alignment by clearing traceback and setting dimensions */
void ALIGNMENT_Reuse(ALIGNMENT*  aln,
                     int         Q,
                     int         T )
{
   aln->N   = 0;
   aln->beg = -1;
   aln->end = -1;
   
   aln->Q   = Q;
   aln->T   = T;
}

/* push trace onto end of alignment */
void ALIGNMENT_Pushback(ALIGNMENT* aln,
                        TRACE*     tr)
{
   /* if debugging, do edgechecks */
   #if DEBUG
      /* if normal state, check bounds of normal dp matrix */
      /* if special state, check bounds of special dp matrix */
   #endif

   aln->traces[aln->N] = *tr;
   aln->N++;

   /* if array is full, resize */
   if (aln->N >= aln->Nalloc - 1) {
      ALIGNMENT_Resize(aln, aln->Nalloc * 2);
   }
}

/* resize TRACE array in ALIGNMENT */
void ALIGNMENT_Resize(ALIGNMENT* aln,
                      int        size)
{
   aln->traces = (TRACE*) realloc(aln->traces, sizeof(TRACE) * size);
   if (aln->traces == NULL) {
      const char* obj_name = "ALIGNMENT";
      fprintf(stderr, "ERROR: Couldn't realloc TRACES for %s: <%p>.\n", obj_name, aln);
      exit(EXIT_FAILURE);
   }
   aln->Nalloc = size;
}

/* Empty ALIGNMENT Array */
void ALIGNMENT_Clear(ALIGNMENT* aln)
{
   aln->N = 0;
}

/* outputs ALIGNMENT to FILE pointer */
void ALIGNMENT_Save(ALIGNMENT* aln,
                    char*      _filename_)
{
   FILE* fp = fopen(_filename_, "w");
   ALIGNMENT_Dump(aln, fp);
   fclose(fp);
}

/* outputs ALIGNMENT to FILE pointer */
void ALIGNMENT_Dump(ALIGNMENT* aln,
                    FILE*      fp)
{
   /* test for bad file pointer */
   if (fp == NULL) {
      const char* obj_name = "ALIGNMENT";
      fprintf(stderr, "ERROR: Bad FILE pointer for printing %s: <%p>.\n", obj_name, aln);
      exit(EXIT_FAILURE);
      return;
   }
   
   fprintf(fp, "# ALIGNMENT (length=%d)\n", aln->N );
   for (unsigned int i = 0; i < aln->N; ++i)
   {
      int st = aln->traces[i].st;
      fprintf(fp, "[%d](%s,%d,%d)\n", i, STATE_NAMES[st], aln->traces[i].i, aln->traces[i].j);
   }
}
