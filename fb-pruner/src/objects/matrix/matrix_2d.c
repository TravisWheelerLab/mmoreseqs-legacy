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
#include "structs.h"
#include "utilities.h"
#include "objects.h"

/* header */
#include "matrix_2d.h"

/* constructor */
MATRIX_2D* MATRIX_2D_Create(  int  R,
                              int  C )
{
   MATRIX_2D* mx = NULL;

   if ( R <= 0 || C <= 0 ) {
      fprintf(stderr, "ERROR: MATRIX_2D Rows and Columns must be a positive size.\n");
      exit(EXIT_FAILURE);
   }

   mx = (MATRIX_2D*) ERRORCHECK_malloc( sizeof(MATRIX_2D), __FILE__, __LINE__, __FUNCTION__ );
   if (mx == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc MATRIX_3D.\n");
      exit(EXIT_FAILURE);
   }

   mx->R       = 0;
   mx->C       = 0;
   mx->Nalloc  = 0;
   mx->data    = NULL;

   MATRIX_2D_Resize( mx, R, C );
   mx->clean = false;

   return mx;
}

/* constructor with all values set to -INF */
MATRIX_2D* MATRIX_2D_Create_Clean(  int  R,
                                    int  C )
{
   MATRIX_2D* mx;
   mx = MATRIX_2D_Create( R, C );
   MATRIX_2D_Clean( mx );
   mx->clean = true;
   return mx;
}

/* destructor */
void* MATRIX_2D_Destroy(MATRIX_2D*  mx)
{
   if (mx == NULL) return mx;

   free(mx->data);
   free(mx);
   mx = NULL;
   return mx;
}

/* deep copy: returns dest matrix, will allocate if null */
MATRIX_2D* MATRIX_2D_Copy( MATRIX_2D*     dest,
                           MATRIX_2D*     src )
{
   /* dim */
   int R = src->R;
   int C = src->C;

   /* create a matrix if necessary */
   if ( dest == NULL ) {
      dest = MATRIX_2D_Create( R, C );
   } else {
      MATRIX_2D_Reuse( dest, R, C );
   }
   /* copy values */
   for ( int i = 0; i < R; i++ ) {
      for ( int j = 0; j < C; j++ ) {
            MX_2D( dest, i, j ) = MX_2D( src, i, j );
      }
   }

   return dest;
}

/* fill MATRIX_2D with values */
void MATRIX_2D_Fill( MATRIX_2D*     mx,
                     float          val)
{
   for (int i = 0; i < mx->R; i++) {
      for (int j = 0; j < mx->C; j++) {
         *MATRIX_2D_Get(mx, i, j) = val;
      }
   }
}

/* fill MATRIX_2D with -INF */
void MATRIX_2D_Clean( MATRIX_2D*   mx)
{
   MATRIX_2D_Fill( mx, -INF );
   mx->clean = true;
}

/* check that all cells are filled with given value */
int MATRIX_2D_Check_Value(  MATRIX_2D*     mx,
                           float          val)
{
   int cnt = 0;
   for (int i = 0; i < mx->R; i++) {
      for (int j = 0; j < mx->C; j++) {
         if ( *MATRIX_2D_Get(mx, i, j) != val )
            cnt++;
      }
   }
   return cnt;
}

/* fill MATRIX_2D with -INF */
int MATRIX_2D_Check_Clean( MATRIX_2D*   mx)
{
   return MATRIX_2D_Check_Value( mx, -INF );
}


/* getter pointer for index in MATRIX */
inline
float* MATRIX_2D_Get(MATRIX_2D*  mx,
                     int         i,
                     int         j)
{
   /* if debugging, do edgechecks */
   #if DEBUG
      int n = MATRIX_2D_to_1D(mx, i, j);
      int used = mx->R * mx->C;
      if (i >= mx->R || i < 0 || j >= mx->C || j < 0 || n >= used ) {
         fprintf(stderr, "ERROR: MATRIX_2D Access Out-of-Bounds\n");
         fprintf(stderr, "2D => dim: (%d,%d), access: (%d,%d)\n", mx->R, mx->C, i, j);
         fprintf(stderr, "1D => dim: (%d/%d), access: (%d)\n", used, mx->Nalloc, n);
         exit(EXIT_FAILURE);
      }
   #endif

   float* data = &( mx->data[ MATRIX_2D_to_1D(mx, i, j) ] );
   return data;
}

/* getter pointer for index in MATRIX (input in 1D-coords) */
inline
float* MATRIX_2D_Get_1D(   MATRIX_2D*  mx,
                           int         n )
{
   // /* if debugging, do edgechecks */
   // #if DEBUG
   //    int used = mx->R * mx->C;
   //    if ( n < 0|| n >= used ) {
   //       fprintf(stderr, "ERROR: MATRIX_2D Access Out-of-Bounds\n");
   //       fprintf(stderr, "1D => dim: (%d/%d), access: (%d)\n", used, mx->Nalloc, n);
   //       exit(EXIT_FAILURE);
   //    }
   // #endif

   float* data = &( mx->data[ n ] );
   return data;
}

/* convert 2D-coords to 1D-coords */
inline
int MATRIX_2D_to_1D(MATRIX_2D*  mx,
                    int         i,
                    int         j)
{
   // assert(i < mx->R && j < mx->C);
   return (i * mx->C + j);
}

/* reuse MATRIX_2D by resizing only if new matrix requires more memory */
float MATRIX_2D_Reuse(  MATRIX_2D*  mx,
                        int         R,
                        int         C )
{
   if (R * C > mx->Nalloc) {
      MATRIX_2D_Resize( mx, R, C );
   }
   else
   {
      mx->R = R;
      mx->C = C;
   }
   mx->clean = false;
}

/* reuse MATRIX_2D by resizing only if new matrix requires more memory.  All new matrix values are cleaned. */
float MATRIX_2D_Reuse_Clean(  MATRIX_2D*  mx,
                              int         R,
                              int         C )
{
   if ( mx->clean == false )
      MATRIX_2D_Clean( mx );

   int N_prv = mx->R * mx->C;

   if (R * C > mx->Nalloc) {
      MATRIX_2D_Resize( mx, R, C );
   }
   else
   {
      mx->R = R;
      mx->C = C;
   }

   for ( int i = N_prv; i < mx->Nalloc; i++ ) {
      *MATRIX_2D_Get_1D( mx, i ) = -INF;
   }
   mx->clean = true;
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
void MATRIX_2D_Dump( MATRIX_2D*  mx,
                     FILE*       fp )
{
   /* check for bad pointer */
   if (fp == NULL) {
      fprintf(stderr, "ERROR: Bad FILE POINTER for printing SEQUENCE.\n");
      exit(EXIT_FAILURE);
   }

   fprintf(fp, "=== MATRIX_2D { R, C } = { %d, %d }===\n", mx->R, mx->C);
   for (int i = 0; i < mx->R; i++) {
      for (int j = 0; j < mx->C; j++) {
         fprintf(fp, "%.1f\t", *MATRIX_2D_Get(mx, i, j) );
      }
      fprintf(fp, "\n");
   }
}

/* Save MATRIX_2D to FILE by FILENAME */
void MATRIX_2D_Save( MATRIX_2D*  mx,
                     char*       _filename_ )
{
   FILE* fp = fopen(_filename_, "w");
   MATRIX_2D_Dump(mx, fp);
   printf("MATRIX_2D saved to '%s'\n", _filename_);
   fclose(fp);
}

/* Compare two MATRIX_2D */
int MATRIX_2D_Compare(  MATRIX_2D*  mx_a,
                        MATRIX_2D*  mx_b )
{
   /* inequality value */
   int      cmp   = 0;
   float    diff  = 0;
   bool     eq    = false;
   /* set float equality tolerance */
   const float tol = 1e-2;

   /* dim */
   int   R  = mx_a->R;
   int   C  = mx_a->C;

   MATRIX_2D* cloud_MX = debugger->cloud_MX;
   #if DEBUG
   {
      MATRIX_2D_Reuse( cloud_MX, R, C );
      MATRIX_2D_Fill( cloud_MX, 0 );
   }
   #endif

   if (  mx_a->R != mx_b->R || 
         mx_a->C != mx_b->C ) {
      return -1;
   }

   for ( int i = 0; i < R; i++ ) {
      for ( int j = 0; j < C; j++ ) {
         eq    = ( MX_2D( mx_a, i, j ) == MX_2D( mx_b, i, j ) );
         diff  = MX_2D( mx_a, i, j ) - MX_2D( mx_b, i, j );
         diff  = ABS( diff );
         if ( eq == false && diff > tol ) {
            #if DEBUG
            {
               if ( cmp == 0 ) {  /* only reports first inequality */
                  printf("MATRIX_2D EQUALITY failed at: (%d,%d), TOLERANCE: %f\n", i, j, tol);
                  printf("MX_A: %f, MX_B: %f => DIFF: %f\n", MX_2D( mx_a, i, j ), MX_2D( mx_b, i, j ), diff );
               }
               MX_2D( cloud_MX, i, j ) += 1.0;
            }
            #endif
            cmp++;
            #if !DEBUG 
            {
               return cmp;
            }
            #endif
         }
      }
   }
   return cmp;
}

/*
 *  FUNCTION:  MATRIX_2D_Diff()
 *  SYNOPSIS:  Takes difference of <mx_a> - <mx_b>.  Result stored in <mx_diff>.
 */
int MATRIX_2D_Diff(  MATRIX_2D*  mx_a,
                     MATRIX_2D*  mx_b,
                     MATRIX_2D*  mx_diff )
{
   for ( int i = 0; i < mx_a->R; i++ ) {
      for ( int j = 0; j < mx_a->C; j++ ) {
         if ( MX_2D( mx_a, i, j ) == MX_2D( mx_b, i, j ) )
            MX_2D( mx_diff, i, j ) = 0;   /* verifies equality when INF or -INF */
         else
            MX_2D( mx_diff, i, j ) = MX_2D( mx_a, i, j ) - MX_2D( mx_b, i, j );
      }
   }
   return 0;
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
         float data = (i + 1) * 100 + (j + 1);
         printf("Inserting %.1f into index (%d,%d)=%d...\n", data, i, j, MATRIX_2D_to_1D(mx, i, j) );
         *MATRIX_2D_Get(mx, i, j) = data;
      }
   }
   MATRIX_2D_Dump(mx, stdout);
}