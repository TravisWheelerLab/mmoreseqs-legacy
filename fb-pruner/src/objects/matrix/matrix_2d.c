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


/*! FUNCTION:  MATRIX_2D_Create()
 *  SYNOPSIS:  Constructor: Allocates memory for MATRIX_2D.
 *
 *  RETURN:    Return pointer to new MATRIX_2D object.
 */
MATRIX_2D* 
MATRIX_2D_Create( int  R,
                  int  C )
{
   MATRIX_2D* mx = NULL;

   if ( R <= 0 || C <= 0 ) {
      fprintf(stderr, "ERROR: MATRIX_2D Rows and Columns must be a positive size.\n");
      exit(EXIT_FAILURE);
   }

   mx = (MATRIX_2D*) ERROR_malloc( sizeof(MATRIX_2D) );

   mx->R       = 0;
   mx->C       = 0;
   mx->Nalloc  = 0;
   mx->data    = NULL;

   MATRIX_2D_Resize( mx, R, C );
   mx->clean = false;

   return mx;
}

/*! FUNCTION:  MATRIX_2D_Create_Clean()
 *  SYNOPSIS:  Constructor: Allocates memory for MATRIX_2D. Sets all data to -INF.
 *
 *  RETURN:    Return pointer to new MATRIX_2D object.
 */
MATRIX_2D* 
MATRIX_2D_Create_Clean(    int  R,
                           int  C )
{
   MATRIX_2D* mx = NULL;

   mx = MATRIX_2D_Create( R, C );
   MATRIX_2D_Clean( mx );
   mx->clean = true;

   return mx;
}

/*! FUNCTION:  MATRIX_2D_Destroy()
 *  SYNOPSIS:  Destructor: Frees memory for MATRIX_2D.
 *
 *  RETURN:    Return pointer to new MATRIX_2D object.
 */
MATRIX_2D* 
MATRIX_2D_Destroy( MATRIX_2D*  mx )
{
   if (mx == NULL) return mx;

   ERROR_free(mx->data);
   ERROR_free(mx);
   mx = NULL;
   return mx;
}

/*! FUNCTION:  MATRIX_2D_Copy()
 *  SYNOPSIS:  Deep Copy: Makes deep copy of <src> to <dest>.  Creates new MATRIX_2D object if <dest> is NULL.
 *
 *  RETURN:    Return pointer to <dest>.
 */
MATRIX_2D* 
MATRIX_2D_Copy( MATRIX_2D*           dest,
                const MATRIX_2D*     src )
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

/*! FUNCTION:  MATRIX_2D_Destroy()
 *  SYNOPSIS:  Destructor: Frees memory for MATRIX_2D.
 *
 *  RETURN:    Return pointer to new MATRIX_2D object.
 */
int 
MATRIX_2D_Fill( MATRIX_2D*  mx,
                float       val )
{
   for (int i = 0; i < mx->R; i++) {
      for (int j = 0; j < mx->C; j++) {
         *MATRIX_2D_Get(mx, i, j) = val;
      }
   }
}

/* fill MATRIX_2D with -INF */
/*! FUNCTION:  
 *  SYNOPSIS:  
 *
 *  RETURN:    
 */
int 
MATRIX_2D_Clean( MATRIX_2D*   mx)
{
   MATRIX_2D_Fill( mx, -INF );
   mx->clean = true;
}

/* check that all cells are filled with given value */
/*! FUNCTION:  
 *  SYNOPSIS:  
 *
 *  RETURN:    
 */
int 
MATRIX_2D_Check_Value(  MATRIX_2D*     mx,
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
/*! FUNCTION:  
 *  SYNOPSIS:  
 *
 *  RETURN:    
 */
int 
MATRIX_2D_Check_Clean( MATRIX_2D*   mx)
{
   return MATRIX_2D_Check_Value( mx, -INF );
}


/* get a pointer to index <i,j> in <mx> */
/*! FUNCTION:  
 *  SYNOPSIS:  
 *
 *  RETURN:    
 */
inline
float* MATRIX_2D_Get(   MATRIX_2D*  mx,
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

/* get a pointer to index <n> in <mx> (input in 1D-coords) */
/*! FUNCTION:  
 *  SYNOPSIS:  
 *
 *  RETURN:    
 */
inline
float* 
MATRIX_2D_Get_1D( MATRIX_2D*  mx,
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
/*! FUNCTION:  
 *  SYNOPSIS:  
 *
 *  RETURN:    
 */
inline
int 
MATRIX_2D_to_1D(  MATRIX_2D*  mx,
                  int         i,
                  int         j)
{
   // assert(i < mx->R && j < mx->C);
   return (i * mx->C + j);
}

/* reuse MATRIX_2D by resizing only if new matrix requires more memory */
/*! FUNCTION:  
 *  SYNOPSIS:  
 *
 *  RETURN:    
 */
int 
MATRIX_2D_Reuse(  MATRIX_2D*  mx,
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
/*! FUNCTION:  
 *  SYNOPSIS:  
 *
 *  RETURN:    
 */
int 
MATRIX_2D_Reuse_Clean(  MATRIX_2D*  mx,
                        int         R,
                        int         C )
{
   /* TODO: fix this function */
   MATRIX_2D_Reuse( mx, R, C );
   MATRIX_2D_Clean( mx );
   return 0.0;

   /* if matrix currently isn't cleaned, do it now */
   if ( mx->clean == false ) {
      MATRIX_2D_Clean( mx );
   }
   /* previous flat array size */
   int N_prv = mx->R * mx->C;

   /* resize to new dimensions */
   if (R * C > mx->Nalloc) {
      MATRIX_2D_Resize( mx, R, C );
   }
   else {
      mx->R = R;
      mx->C = C;
   }
   /* all cells that have been added to matrix must be cleaned */
   for ( int i = N_prv; i < mx->Nalloc; i++ ) {
      *MATRIX_2D_Get_1D( mx, i ) = -INF;
   }
   mx->clean = true;

   return STATUS_SUCCESS;
}


/* resize MATRIX_2D to new dimensions */
/*! FUNCTION:  
 *  SYNOPSIS:  
 *
 *  RETURN:    
 */
int 
MATRIX_2D_Resize( MATRIX_2D*  mx,
                  int         R,
                  int         C)
{
   mx->Nalloc  = R * C;
   mx->R       = R;
   mx->C       = C;

   mx->data = (float*) ERROR_realloc( mx->data, sizeof(float) * (R * C) );
}

/* Outputs MATRIX_2D out to FILE POINTER */
/*! FUNCTION:  
 *  SYNOPSIS:  
 *
 *  RETURN:    
 */
int 
MATRIX_2D_Dump(   MATRIX_2D*  mx,
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
/*! FUNCTION:  
 *  SYNOPSIS:  
 *
 *  RETURN:    
 */
int 
MATRIX_2D_Save(   MATRIX_2D*  mx,
                  char*       _filename_ )
{
   FILE* fp = fopen(_filename_, "w");
   MATRIX_2D_Dump(mx, fp);
   printf("MATRIX_2D saved to '%s'\n", _filename_);
   fclose(fp);
}

/* Compare two MATRIX_2D */
/*! FUNCTION:  
 *  SYNOPSIS:  
 *
 *  RETURN:    
 */
int 
MATRIX_2D_Compare(   MATRIX_2D*  mx_A,
                     MATRIX_2D*  mx_B )
{
   /* inequality value */
   int      cmp   = 0;
   float    diff  = 0;
   bool     eq    = false;
   /* set float equality tolerance */
   const float tol = 1e-2;

   /* dim */
   int   R  = mx_A->R;
   int   C  = mx_A->C;

   MATRIX_2D* cloud_MX = debugger->cloud_MX;
   #if DEBUG
   {
      MATRIX_2D_Reuse( cloud_MX, R, C );
      MATRIX_2D_Fill( cloud_MX, 0 );
   }
   #endif

   if (  mx_A->R != mx_B->R || 
         mx_A->C != mx_B->C ) {
      return -1;
   }

   for ( int i = 0; i < R; i++ ) {
      for ( int j = 0; j < C; j++ ) {
         eq    = ( MX_2D( mx_A, i, j ) == MX_2D( mx_B, i, j ) );
         diff  = MX_2D( mx_A, i, j ) - MX_2D( mx_B, i, j );
         diff  = ABS( diff );
         if ( eq == false && diff > tol ) {
            #if DEBUG
            {
               if ( cmp == 0 ) {  /* only reports first inequality */
                  printf("MATRIX_2D EQUALITY failed at: (%d,%d), TOLERANCE: %f\n", i, j, tol);
                  printf("MX_A: %f, MX_B: %f => DIFF: %f\n", MX_2D( mx_A, i, j ), MX_2D( mx_B, i, j ), diff );
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

/*! FUNCTION:  
 *  SYNOPSIS:  
 *
 *  RETURN:    
 */
int MATRIX_2D_Add(   MATRIX_2D*  mx_A,
                     MATRIX_2D*  mx_B,
                     MATRIX_2D*  mx_res )
{
   for ( int i = 0; i < mx_A->R; i++ ) {
      for ( int j = 0; j < mx_A->C; j++ ) {
         if ( MX_2D( mx_A, i, j ) == -INF || MX_2D( mx_B, i, j ) == -INF ) {
            MX_2D( mx_res, i, j ) = -INF;
         }
         else {
            MX_2D( mx_res, i, j ) = MX_2D( mx_A, i, j ) - MX_2D( mx_B, i, j );
         }
      }
   }
   return STATUS_SUCCESS;
}

/*! FUNCTION:  
 *  SYNOPSIS:  
 *
 *  RETURN:    
 */
int 
MATRIX_2D_Diff(   MATRIX_2D*  mx_A,
                  MATRIX_2D*  mx_B,
                  MATRIX_2D*  mx_diff )
{
   for ( int i = 0; i < mx_A->R; i++ ) {
      for ( int j = 0; j < mx_A->C; j++ ) {
         if ( MX_2D( mx_A, i, j ) == MX_2D( mx_B, i, j ) )
            MX_2D( mx_diff, i, j ) = 0;   /* verifies equality when INF or -INF */
         else
            MX_2D( mx_diff, i, j ) = MX_2D( mx_A, i, j ) - MX_2D( mx_B, i, j );
      }
   }
   return 0;
}

/*! FUNCTION:  
 *  SYNOPSIS:  
 *
 *  RETURN:    
 */
int 
MATRIX_2D_UnitTest()
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