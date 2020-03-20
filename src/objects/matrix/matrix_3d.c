/*******************************************************************************
 *  @file matrix.c
 *  @brief MATRIX_3D Float object.
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
#include <ctype.h>
#include <time.h>

/* local imports */
#include "../structs.h"
#include "../../utility.h"

/* header */
#include "matrix_3d.h"

/* constructor */
MATRIX_3D* MATRIX_3D_Create(int  R, 
                            int  C,
                            int  N) 
{
   MATRIX_3D* mx;

   mx = (MATRIX_3D*) malloc( sizeof(MATRIX_3D) );
   if (mx == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc MATRIX_3D.\n");
      exit(EXIT_FAILURE);
   }

   MATRIX_3D_Resize( mx, R, C, N );

   return mx;
}

/* destructor */
void MATRIX_3D_Destroy(MATRIX_3D*  mx)
{
   free(mx->data);
   free(mx);
}

/* fill MATRIX_3D with values */
MATRIX_3D* MATRIX_3D_Fill(MATRIX_3D*  mx,
                          int         val)
{
   for (int i = 0; i < mx->R; i++) {
      for (int j = 0; j < mx->C; j++) {
         for (int k = 0; k < mx->N; k++) {
            *MATRIX_3D_Get(mx, i, j, k) = val;
         }
      }
   }
}

/* getter for index */
inline
float* MATRIX_3D_Get(MATRIX_3D*  mx, 
                     int         i, 
                     int         j,
                     int         k)
{
   float* data = &( mx->data[ MATRIX_3D_to_1D(mx, i, j, k) ] );
   return data;
}

/* convert 3D-coords to 1D-coords */
inline
int MATRIX_3D_to_1D(MATRIX_3D*  mx, 
                    int         i, 
                    int         j,
                    int         k)
{
   // assert(i < mx->R && j < mx->C);
   return (i * (mx->C * mx->R) + j * (mx->C) + k);
}

/* reuse MATRIX_3D by resizing only if new matrix requires more memory */
float MATRIX_3D_Reuse(MATRIX_3D*  mx, 
                      int         R, 
                      int         C,
                      int         N)
{
   if (R*C*N > mx->Nalloc) {
      MATRIX_3D_Resize(mx, R, C, N);
   } 
   else
   {
      mx->R = R;
      mx->C = C;
      mx->N = N;
   }
}

/* resize MATRIX_3D to new dimensions */
float MATRIX_3D_Resize(MATRIX_3D*  mx, 
                       int         R, 
                       int         C,
                       int         N)
{
   mx->Nalloc = R * C * N;
   mx->data = (float*) realloc( mx->data, sizeof(float) * (R * C * N) );
   if (mx->data == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc DATA for MATRIX_3D.\n");
      exit(EXIT_FAILURE);
   }
   mx->R = R;
   mx->C = C;
   mx->N = N;
}
/* Outputs MATRIX_3D out to FILE POINTER */
void MATRIX_3D_Dump(MATRIX_3D*  mx,
                    FILE*       fp)
{
   /* check for bad pointer */
   if (fp == NULL) {
      fprintf(stderr, "ERROR: Bad FILE POINTER for printing SEQUENCE.\n");
      exit(EXIT_FAILURE);
   }

   fprintf(fp, "=== MATRIX_3D ===\n");
   for (int i = 0; i < mx->N; i++) {
      for (int j = 0; j < mx->R; j++) {
         for (int k = 0; k < mx->C; k++) {
            fprintf(fp, "%.1f\t", *MATRIX_3D_Get(mx, i, j, k) );
         }
         fprintf(fp, "\n");
      }
      fprintf(fp, "\n");
   }
}

/* Save MATRIX_3D to FILE by FILENAME */
void MATRIX_3D_Save(MATRIX_3D*  mx,
                    char*       _filename_)
{
   FILE* fp = fopen(_filename_, "w");
   MATRIX_3D_Dump(mx, fp);
   fclose(fp);
}

/* unit test */
void MATRIX_3D_UnitTest()
{
   int R = 10;
   int C = 5;
   int N = 3;
   MATRIX_3D* mx = MATRIX_3D_Create(R, C, N);

   for (int i = 0; i < R; i++) {
      for (int j = 0; j < C; j++) {
         for (int k = 0; k < N; k++) {
            *MATRIX_3D_Get(mx, i, j, k) = (float)(i*10000 + j*100 + k);
         }
      }
   }
   MATRIX_3D_Dump(mx, stdout);
}