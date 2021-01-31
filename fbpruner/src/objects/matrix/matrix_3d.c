/*******************************************************************************
 *  FILE:      matrix_3d.c
 *  PURPOSE:   MATRIX_3D Float object.
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
#include "../../utilities/_utilities.h"
#include "../_objects.h"

/* header */
#include "matrix_3d.h"

/* constructor */
MATRIX_3D* 
MATRIX_3D_Create( const int  R,
                  const int  C,
                  const int  N )
{
   if ( R <= 0 || C <= 0 ) {
      fprintf(stderr, "ERROR: MATRIX_2D Rows and Columns must be a positive size.\n");
      exit(EXIT_FAILURE);
   }

   MATRIX_3D* mx;
   mx = ERROR_malloc( sizeof(MATRIX_3D) );

   mx->R       = 0;
   mx->C       = 0;
   mx->N       = 0;
   mx->Nalloc  = 0;
   mx->data    = NULL;
   mx->clean   = false;

   MATRIX_3D_Reuse( mx, R, C, N );

   return mx;
}

/* constructor that allocates no matrix data, leaves NULL to be supplied by user */
MATRIX_3D* 
MATRIX_3D_Create_NoData(   const int  R,
                           const int  C,
                           const int  N )
{
   if ( R <= 0 || C <= 0 ) {
      fprintf(stderr, "ERROR: MATRIX_2D Rows and Columns must be a positive size.\n");
      exit(EXIT_FAILURE);
   }

   MATRIX_3D* mx;
   mx = ERROR_malloc( sizeof(MATRIX_3D) );

   mx->R       = 0;
   mx->C       = 0;
   mx->N       = 0;
   mx->Nalloc  = 0;
   mx->data    = NULL;
   mx->clean   = false;

   MATRIX_3D_Reuse( mx, R, C, N );

   return mx;
}

/* constructor for clean matrices */
MATRIX_3D* 
MATRIX_3D_Create_Clean( const int  R,
                        const int  C,
                        const int  N )
{
   MATRIX_3D* mx;

   mx = MATRIX_3D_Create( R, C, N );
   MATRIX_3D_Clean( mx );
   mx->clean = true;
   return mx;
}

/* destructor */
MATRIX_3D* 
MATRIX_3D_Destroy( MATRIX_3D*  mx )
{
   if (mx == NULL) return NULL;

   ERROR_free(mx->data);
   ERROR_free(mx);

   return NULL;
}

/* deep copy: returns dest matrix, will allocate if null */
MATRIX_3D* 
MATRIX_3D_Copy(   MATRIX_3D*           dest,
                  MATRIX_3D*           src )
{
   /* verify dest and src have same dimensions */
   #if DEBUG
   {
      if ( dest->R != src->R || dest->C != src->C || dest->N != src->N ) {
         printf("ERROR: src and dest do not have the same dimensions.\n");
         exit(EXIT_FAILURE);
      }
   }
   #endif

   /* dim */
   int R = src->R;
   int C = src->C;
   int N = src->N;

   /* create a matrix if necessary */
   if ( dest == NULL ) {
      dest = MATRIX_3D_Create( R, C, N );
   } else {
      MATRIX_3D_Reuse( dest, R, C, N );
   }
   /* copy values */
   for ( int i = 0; i < R; i++ ) {
      for ( int j = 0; j < C; j++ ) {
         for ( int k = 0; k < N; k++ ) {
            *MATRIX_3D_Get( dest, i, j, k ) = *MATRIX_3D_Get( src, i, j, k );
         }
      }
   }

   return dest;
}

/* fill MATRIX_3D with values */
void MATRIX_3D_Fill( MATRIX_3D*     mx,
                     const float    val )
{
   for (int i = 0; i < mx->R; i++) {
      for (int j = 0; j < mx->C; j++) {
         for (int k = 0; k < mx->N; k++) {
            *MATRIX_3D_Get(mx, i, j, k) = val;
         }
      }
   }
}


/* fill MATRIX_3D with -INF */
void MATRIX_3D_Clean(   MATRIX_3D*     mx )
{
   MATRIX_3D_Fill( mx, -INF );
   mx->clean = true;
}

/* check MATRIX_3D and counts all cells not matching val */
int MATRIX_3D_Check_Value( MATRIX_3D*     mx,
                           const float    val )
{
   int cnt = 0;
   int N = mx->R * mx->C * mx->N;
   for ( int i = 0; i < mx->R; i++ )
      for ( int j = 0; j < mx->C; j++ )
         for ( int k = 0; k < mx->N; k++ )
            if ( *MATRIX_3D_Get( mx, i, j, k ) != val ) {
               cnt++; 
               /* NOTE: turn on? */
               // printf("#> CHECK VALUE ERROR: at (%d,%d,%d) = %9.4f \n", i,j,k, *MATRIX_3D_Get(mx, i,j,k) );
            }

   return cnt;
}

/* check MATRIX_3D and verify all cells contain -INF */
int MATRIX_3D_Check_Clean( MATRIX_3D*     mx )
{
   return MATRIX_3D_Check_Value( mx, -INF );
}

/* getter for index */
inline
float* MATRIX_3D_Get(   MATRIX_3D*  mx,
                        const int   i,
                        const int   j,
                        const int   k )
{
   /* if debugging, do edgebound checks */
   #if DEBUG
   {
      int n = MATRIX_3D_to_1D(mx, i, j, k);
      int used = mx->R * mx->C * mx->N;
      if (i >= mx->R || i < 0 || j >= mx->C || j < 0 || k >= mx->N || n >= used ) {
         fprintf(stderr, "ERROR: MATRIX_3D Access Out-of-Bounds\n");
         fprintf(stderr, "3D => dim: (%d,%d,%d), access: (%d,%d,%d)\n", mx->R, mx->C, mx->N, i, j, k);
         fprintf(stderr, "1D => dim: (%d/%d), access: (%d)\n", used, mx->Nalloc, n);
         exit(EXIT_FAILURE);
      }
   }
   #endif

   float* data = &( mx->data[ MATRIX_3D_to_1D(mx, i, j, k) ] );
   return data;
}

/** FUNCTION:  MATRIX_3D_Get_X()
 *  SYNOPSIS:  Get pointer to cell from <mx> at position <i,j,k>.
 */
inline
float* 
MATRIX_3D_Get_X(  MATRIX_3D*  mx,
                  const int   i,
                  const int   j,
                  const int   k )
{
   /* if debugging, do edgebound checks */
   #if DEBUG
   int n = MATRIX_3D_to_1D(mx, i, j, k);
   int used = mx->R * mx->C * mx->N;
   if (i >= mx->R || i < 0 || j >= mx->C || j < 0 || k >= mx->N || n >= used ) {
      fprintf(stderr, "ERROR: MATRIX_3D Access Out-of-Bounds\n");
      fprintf(stderr, "3D => dim: (%d,%d,%d), access: (%d,%d,%d)\n", mx->R, mx->C, mx->N, i, j, k);
      fprintf(stderr, "1D => dim: (%d/%d), access: (%d)\n", used, mx->Nalloc, n);
      exit(EXIT_FAILURE);
   }
   #endif

   float* data = &( mx->data[ MATRIX_3D_to_1D(mx, i, j, k) ] );
   return data;
}


/** FUNCTION:  MATRIX_3D_Get_1D()
 *  SYNOPSIS:  Get pointer to cell from <mx> at flat-index position <n>.
 */
float* 
MATRIX_3D_Get_1D( MATRIX_3D*  mx,
                  const int   n )
{
   /* if debugging, do edgebound checks */
   #if DEBUG
   int used = mx->R * mx->C * mx->N;
   if ( n < 0 || n >= used ) {
      fprintf(stderr, "ERROR: MATRIX_3D Access Out-of-Bounds\n");
      fprintf(stderr, "1D => dim: (%d/%d), access: (%d)\n", used, mx->Nalloc, n);
      exit(EXIT_FAILURE);
   }
   #endif

   float* data = &( mx->data[n] );
   return data;
}


/** FUNCTION:  MATRIX_3D_to_1D()
 *  SYNOPSIS:  
 */
inline
int 
MATRIX_3D_to_1D(  const MATRIX_3D*  mx,
                  const int         i,
                  const int         j,
                  const int         k )
{
   /* (i,j,k) -> (R,C,N) */
   return ((i * (mx->C * mx->N)) + (j * (mx->N)) + k);
}


/** FUNCTION:  MATRIX_3D_Reuse()
 *  SYNOPSIS:  Reuse MATRIX_3D by resizing only if new matrix requires more memory.
 */
float 
MATRIX_3D_Reuse(  MATRIX_3D*  mx,
                  const int   R,
                  const int   C,
                  const int   N )
{
   if (R * C * N > mx->Nalloc) {
      MATRIX_3D_Resize(mx, R, C, N);
   }
   else
   {
      mx->R = R;
      mx->C = C;
      mx->N = N;
   }
   mx->clean = false;
}


/** FUNCTION:  MATRIX_3D_Reuse_Clean()
 *  SYNOPSIS:  Reuse MATRIX_3D by resizing only if new matrix requires more memory. 
 *             New memory is set to -INF.
 */
float 
MATRIX_3D_Reuse_Clean(  MATRIX_3D*  mx,
                        const int   R,
                        const int   C,
                        const int   N )
{
   /* TODO: fix this function */
   MATRIX_3D_Reuse( mx, R, C, N );
   MATRIX_3D_Clean( mx );
   return 0.0;

   if ( mx->clean == false ) {
      MATRIX_3D_Clean( mx );
   }
   int prv_dim = mx->R * mx->C * mx->N;
   int new_dim = R * C * N;

   MATRIX_3D_Reuse( mx, R, C, N );

   for ( int i = prv_dim; i < new_dim; i++ ) {
      *MATRIX_3D_Get_1D( mx, i ) = -INF;
   }
   mx->clean = true;

   #if DEBUG
   {
      MATRIX_3D_Clean( mx );
   }
   #endif
}


/** FUNCTION:  MATRIX_3D_Resize()
 *  SYNOPSIS:  Resize MATRIX_3D to new dimensions.
 */
float 
MATRIX_3D_Resize( MATRIX_3D*  mx,
                  const int   R,
                  const int   C,
                  const int   N )
{
   mx->Nalloc  = R * C * N;
   mx->data    = ERROR_realloc( mx->data, sizeof(float) * (R * C * N) );
   mx->R       = R;
   mx->C       = C;
   mx->N       = N;
}


/** FUNCTION:  MATRIX_3D_Dump()
 *  SYNOPSIS:  Outputs MATRIX_3D <mx> out to file pointer <fp>.
 */
void 
MATRIX_3D_Dump(   MATRIX_3D*  mx,
                  FILE*       fp )
{
   /* check for bad pointer */
   if (fp == NULL) {
      fprintf(stderr, "ERROR: Bad FILE POINTER for printing SEQUENCE.\n");
      exit(EXIT_FAILURE);
   }

   fprintf(fp, "=== MATRIX_3D { R, C, N } { %d, %d, %d } ===\n", mx->R, mx->C, mx->N);
   for (int i = 0; i < mx->R; i++) {
      for (int j = 0; j < mx->C; j++) {
         for (int k = 0; k < mx->N; k++) {
            fprintf(fp, "%8.4f\t", *MATRIX_3D_Get(mx, i, j, k) );
         }
         fprintf(fp, "\n");
      }
      fprintf(fp, "\n");
   }
}


/** FUNCTION:  MATRIX_3D_Save()
 *  SYNOPSIS:  Save MATRIX_3D <mx> to file with filename <filename>.
 */
void 
MATRIX_3D_Save(   MATRIX_3D*  mx,
                  char*       _filename_)
{
   FILE* fp = fopen(_filename_, "w");
   MATRIX_3D_Dump(mx, fp);
   fclose(fp);
}


/** FUNCTION:  MATRIX_3D_Compare()
 *  SYNOPSIS:  Compare MATRIX_3D's <mx_A> and <mx_B>.
 */
int 
MATRIX_3D_Compare(   MATRIX_3D*  mx_a,
                     MATRIX_3D*  mx_b )
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
   int   N  = mx_a->N;

   MATRIX_2D* cloud_MX = debugger->cloud_MX;
   #if DEBUG
   {
      MATRIX_2D_Reuse( cloud_MX, C, N );
      MATRIX_2D_Fill( cloud_MX, 0 );
   }
   #endif

   if (  mx_a->R != mx_b->R ||
         mx_a->C != mx_b->C ||
         mx_a->N != mx_b->N  ) {
      return -1;
   }

   for ( int i = 0; i < R; i++ ) {
      for ( int j = 0; j < C; j++ ) {
         for ( int k = 0; k < N; k++ ) {
            /* eq checks for cases INF or -INF */
            eq    = ( MX_3D( mx_a, i, j, k ) == MX_3D( mx_b, i, j, k ) );
            diff  = MX_3D( mx_a, i, j, k ) - MX_3D( mx_b, i, j, k );
            diff  = ABS( diff );
            if ( eq == false && diff > tol ) 
            {
               #if DEBUG
               {
                  if (cmp == 0 ) { /* only reports first inequality */
                     printf("MATRIX_3D EQUALITY failed at: (%d,%d,%d), TOLERANCE: %f\n", i, j, k, tol);
                     printf("MX_A: %f, MX_B: %f => DIFF: %f\n", MX_3D( mx_a, i, j, k ), MX_3D( mx_b, i, j, k ), diff );
                  }
                  MX_2D( cloud_MX, j, k ) += 1.0;
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
   }
   return cmp;
}


/** FUNCTION:  MATRIX_3D_Add()
 *  SYNOPSIS:  Takes sum of <mx_a> + <mx_b>.  Result stored in <mx_res>.
 */
void 
MATRIX_3D_Add( MATRIX_3D*  mx_a,
               MATRIX_3D*  mx_b,
               MATRIX_3D*  mx_res )
{
   for ( int i = 0; i < mx_a->R; i++ ) {
      for ( int j = 0; j < mx_a->C; j++ ) {
         for ( int k = 0; k < mx_a->N; k++ ) {
            MX_3D( mx_res, i, j, k ) = MX_3D( mx_a, i, j, k ) + MX_3D( mx_b, i, j, k );
         }
      }
   }
}


/** FUNCTION:  MATRIX_3D_Diff()
 *  SYNOPSIS:  Takes difference of <mx_a> - <mx_b>.  Result stored in <mx_diff>.
 */
void 
MATRIX_3D_Diff(   MATRIX_3D*  mx_a,
                  MATRIX_3D*  mx_b,
                  MATRIX_3D*  mx_diff )
{
   for ( int i = 0; i < mx_a->R; i++ ) {
      for ( int j = 0; j < mx_a->C; j++ ) {
         for ( int k = 0; k < mx_a->N; k++ ) {
            if ( MX_3D( mx_a, i, j, k ) == MX_3D( mx_b, i, j, k ) )
               MX_3D( mx_diff, i, j, k ) = 0;   /* verifies equality when INF or -INF */
            else
               MX_3D( mx_diff, i, j, k ) = MX_3D( mx_a, i, j, k ) - MX_3D( mx_b, i, j, k );
         }
      }
   }
}


/** FUNCTION:  MATRIX_3D_Log()
 *  SYNOPSIS:  Takes the logrithm function log() of each value in the matrix.
 */
void 
MATRIX_3D_Log(   MATRIX_3D*  mx )
{
   for ( int i = 0; i < mx->R; i++ ) {
      for ( int j = 0; j < mx->C; j++ ) {
         for ( int k = 0; k < mx->N; k++ ) {
            MX_3D( mx, i, j, k ) = log( MX_3D( mx, i, j, k ) );
         }
      }
   }
}


/** FUNCTION:  MATRIX_3D_Exp()
 *  SYNOPSIS:  Takes the exponential function exp() of each value in the matrix.
 */
void 
MATRIX_3D_Exp(   MATRIX_3D*  mx )
{
   for ( int i = 0; i < mx->R; i++ ) {
      for ( int j = 0; j < mx->C; j++ ) {
         for ( int k = 0; k < mx->N; k++ ) {
            MX_3D( mx, i, j, k ) = log( MX_3D( mx, i, j, k ) );
         }
      }
   }
}


/** FUNCTION:  MATRIX_3D_Utest()
 *  SYNOPSIS:  Unit Test for MATIX_3D.
 */
void MATRIX_3D_Utest()
{
   int R = 10;
   int C = 5;
   int N = 3;
   MATRIX_3D* mx = MATRIX_3D_Create(R, C, N);

   for (int i = 0; i < R; i++) {
      for (int j = 0; j < C; j++) {
         for (int k = 0; k < N; k++) {
            *MATRIX_3D_Get(mx, i, j, k) = (float)(i * 10000 + j * 100 + k);
         }
      }
   }
   MATRIX_3D_Dump(mx, stdout);

   MATRIX_3D_Fill(mx, -INF);
   MATRIX_3D_Dump(mx, stdout);

   MATRIX_3D_Destroy(mx);
}