/*******************************************************************************
 *  FILE:    alignment.c
 *  PURPOSE: ALIGNMENT Object.
 *
 *  AUTHOR:  Dave Rich
 *  BUG:     Lots.
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
#include "../../utility.h"

/* header */
#include "alignment.h"

/* constructor */
ALIGNMENT* ALIGNMENT_Create()
{
   const int min_size = 16;
   ALIGNMENT *aln;
   aln = (ALIGNMENT*) malloc( sizeof(ALIGNMENT) );
   if (aln == NULL) {
      fprintf(stderr, "ERROR: Couldn't malloc ALIGNMENT.\n");
      exit(EXIT_FAILURE);
   }
   aln->traces = (TRACE*) malloc( sizeof(TRACE) * min_size );
      if (aln->traces == NULL) {
      fprintf(stderr, "ERROR: Couldn't malloc TRACES for ALIGNMENT.\n");
      exit(EXIT_FAILURE);
   }
}

/* destructor */
void ALIGNMENT_Destroy(ALIGNMENT* aln)
{
   free(aln->traces);
   free(aln);
}

/* push trace onto end of alignment */
void ALIGNMENT_Pushback(ALIGNMENT* aln,
                        TRACE*     tr)
{
   aln->traces[aln->N] = *tr;
   aln->N++;

   /* if array is full, resize */
   if (aln->N >= aln->Nalloc - 1) {
      ALIGNMENT_Resize(aln, 2);
   }
}

/* resize TRACE array in ALIGNMENT */
void ALIGNMENT_Resize(ALIGNMENT* aln,
                      float      growth_factor)
{
   aln->traces = (TRACE*) realloc(aln->traces, sizeof(TRACE) * aln->Nalloc * growth_factor);
   if (aln->traces == NULL) {
      fprintf(stderr, "ERROR: Couldn't realloc TRACES for ALIGNMENT.\n");
      exit(EXIT_FAILURE);
   }
   aln->Nalloc *= growth_factor;
}


/* outputs ALIGNMENT to FILE pointer */
void ALIGNMENT_Dump(ALIGNMENT* aln,
                    FILE*      fp)
{
   /* test for bad file pointer */
   if (fp == NULL) {
      const char* obj_name = "ALIGNMENT";
      fprintf(stderr, "ERROR: Bad FILE pointer for printing %s.\n", obj_name);
      exit(EXIT_FAILURE);
      return;
   }

   static char * states[] = {"ST_M",
                             "ST_I",
                             "ST_D",
                             "ST_E",
                             "ST_N",
                             "ST_J",
                             "ST_C",
                             "ST_B",
                             "ST_S",
                             "ST_T",
                             "ST_X" };
   
   for (unsigned int i = 0; i < aln->N; ++i)
   {
      int st = aln->traces[i].st;
      fprintf(fp, "[%d](%s,%d,%d)\n", i, states[st], aln->traces[i].i, aln->traces[i].j);
   }
}


/* outputs ALIGNMENT to FILE pointer */
void ALIGNMENT_Save(ALIGNMENT* aln,
                    char*      _filename_)
{
   FILE* fp = fopen(_filename_, "w");

   ALIGNMENT_Dump(aln, fp);

   fclose(fp);
}