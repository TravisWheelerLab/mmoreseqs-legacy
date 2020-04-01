/*******************************************************************************
 *  FILE:      matrix_2d.c
 *  PURPOSE:   MATRIX_2D Float object.
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
#include <ctype.h>
#include <time.h>

/* local imports */
#include "../structs.h"
#include "../../utility.h"

/* header */
#include "matrix_2d.h"

/* constructor */
MATRIX_2D* MATRIX_2D_Create(int  R, 
                            int  C) 
{
   MATRIX_2D* mx = NULL;

   if ( R <= 0 || C <= 0 ) {
      fprintf(stderr, "ERROR: MATRIX_2D Rows and Columns must be a positive size.\n");
      exit(EXIT_FAILURE);
   }

   mx = (MATRIX_2D*) malloc( sizeof(MATRIX_2D) );
   if (mx == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc MATRIX_3D.\n");
      exit(EXIT_FAILURE);
   }

   mx->R       = 0;
   mx->C       = 0;
   mx->Nalloc  = 0;
   mx->data    = NULL;

   MATRIX_2D_Resize( mx, R, C );
   
   return mx;
}

/* destructor */
void MATRIX_2D_Destroy(MATRIX_2D*  mx)
{
   free(mx->data);
   free(mx);
}

/* fill MATRIX_2D with values */
MATRIX_2D* MATRIX_2D_Fill(MATRIX_2D*  mx,
                          int         R, 
                          int         C, 
                          int         val)
{
   for (int i = 0; i < R; i++) {
      for (int j = 0; j < C; j++) {
         *MATRIX_2D_Get(mx, i, j) = val;
      }
   }
}

/* getter pointer for index in MATRIX */
inline
float* MATRIX_2D_Get(MATRIX_2D*  mx, 
                     int         i, 
                     int         j)
{
   float* data = &( mx->data[ MATRIX_2D_to_1D(mx, i, j) ] );
   return data;
}

/* convert 2D-coords to 1D-coords */
inline
int MATRIX_2D_to_1D(MATRIX_2D*  mx, 
                    int         i, 
                    int         j)
{
   // assert(i < mx->R && j < mx->C);
   return (i*mx->C + j);
}

/* reuse MATRIX_2D by resizing only if new matrix requires more memory */
float MATRIX_2D_Reuse(MATRIX_2D*  mx, 
                      int         R, 
                      int         C)
{
   if (R*C > mx->Nalloc) {
      MATRIX_2D_Resize(mx, R, C);
   } 
   else
   {
      mx->R = R;
      mx->C = C;
   }
}

/* resize MATRIX_2D to new dimensions */
float MATRIX_2D_Resize(MATRIX_2D*  mx, 
                       int         R, 
                       int         C)
{
   mx->Nalloc  = R * C;
   mx->R       = R;
   mx->C       = C;

   mx->data = (float*) realloc( mx->data, sizeof(float) * (R * C) );
   if (mx->data == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc DATA for MATRIX_3D.\n");
      exit(EXIT_FAILURE);
   }
}
/* Outputs MATRIX_2D out to FILE POINTER */
void MATRIX_2D_Dump(MATRIX_2D*  mx,
                    FILE*       fp)
{
   /* check for bad pointer */
   if (fp == NULL) {
      fprintf(stderr, "ERROR: Bad FILE POINTER for printing SEQUENCE.\n");
      exit(EXIT_FAILURE);
   }

   fprintf(fp, "=== MATRIX_2D ===\n");
   for (int i = 0; i < mx->R; i++) {
      for (int j = 0; j < mx->C; j++) {
         fprintf(fp, "%.1f\t", *MATRIX_2D_Get(mx, i, j) );
      }
      fprintf(fp, "\n");
   }
}

/* Save MATRIX_2D to FILE by FILENAME */
void MATRIX_2D_Save(MATRIX_2D*  mx,
                    char*       _filename_)
{
   FILE* fp = fopen(_filename_, "w");
   MATRIX_2D_Dump(mx, fp);
   fclose(fp);
}

/* unit test */
void MATRIX_2D_UnitTest()
{
   int R = 10;
   int C = 5;
   MATRIX_2D* mx = MATRIX_2D_Create(R, C);
   printf("matrix created...\n");
   printf("MX=>%p\n", mx);

   for (int i = 0; i < R; i++) {
      for (int j = 0; j < C; j++) {
         float data = (i+1)*100 + (j+1);
         printf("Inserting %.1f into index (%d,%d)=%d...\n", data, i, j, MATRIX_2D_to_1D(mx, i, j) );
         *MATRIX_2D_Get(mx, i, j) = data;
      }
   }
   MATRIX_2D_Dump(mx, stdout);
}