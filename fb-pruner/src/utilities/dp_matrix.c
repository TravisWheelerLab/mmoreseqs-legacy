/*******************************************************************************
 *  FILE:      dp_matrix.h
 *  SYNOPSIS:  Dynamic Programming Matrix functions.
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
#include "objects.h"

/* header */
#include "utilities.h"

/*
 *  FUNCTION:  DP_MATRIX_Get_Bounds()
 *  SYNOPSIS:  Get the edgebounds of matrix at given antidiagonal (closed form).
 *             Stores new bounds in BOUND <bnd>.
 */
inline
void DP_MATRIX_Get_Bounds( const int   Q,
                           const int   T,
                           int         d_0,
                           int         d_cnt,
                           BOUND*      bnd )
{

}

/*
 *  FUNCTION:  DP_MATRIX_Copy()
 *  SYNOPSIS:  Copy dynamic programming matrix into destination.
 */
void DP_MATRIX_Copy( const int   Q,
                     const int   T,
                     MATRIX_3D*  st_MX_src,
                     MATRIX_2D*  sp_MX_src,
                     MATRIX_3D*  st_MX_dst,
                     MATRIX_2D*  sp_MX_dst )
{
   int i, j, st;

   for (i = 0; i <= Q; i++)  {
      for (j = 0; j <= T; j++) {
         for (st = 0; st < NUM_NORMAL_STATES; st++) {
            MX_3D(st_MX_dst, st, i, j) = MX_3D(st_MX_src, st, i, j);
         }
      }

      for (st = 0; st < NUM_SPECIAL_STATES; st++) {
         MX_2D(sp_MX_dst, st, i) = MX_2D(sp_MX_src, st, i);
      }
   }
}

/*
 *  FUNCTION:  DP_MATRIX_Fill()
 *  SYNOPSIS:  Fill entire dynamic programming matrix with value
 */
void DP_MATRIX_Fill( const int   Q, 
                     const int   T,
                     MATRIX_3D*  st_MX,
                     MATRIX_2D*  sp_MX,
                     const float val )
{
   MATRIX_3D_Fill( st_MX, val );
   MATRIX_2D_Fill( sp_MX, val );
}

/*
 *  FUNCTION:  DP_MATRIX_Clean()
 *  SYNOPSIS:  Fill entire dynamic programming matrix with clean.
 */
void DP_MATRIX_Clean(   const int   Q, 
                        const int   T,
                        MATRIX_3D*  st_MX,
                        MATRIX_2D*  sp_MX )
{
   MATRIX_3D_Clean( st_MX );
   MATRIX_2D_Clean( sp_MX );
}

/*
 *  FUNCTION:  DP_MATRIX_Clean_Verify()
 *  SYNOPSIS:  Check whether there are clean.  If clean, returns true.
 */
bool DP_MATRIX_Clean_Verify(  const int   Q, 
                              const int   T,
                              MATRIX_3D*  st_MX,
                              MATRIX_2D*  sp_MX )
{
   return st_MX->clean && sp_MX->clean;
}

/*
 *  FUNCTION:  DP_MATRIX_Save()
 *  SYNOPSIS:  Compare two dynamic programming matrices. 
 *    RETURN:  Returns 0 if equal; otherwise returns count of differing cells.
 */
int DP_MATRIX_Compare ( MATRIX_3D*  st_MX_1,
                        MATRIX_2D*  sp_MX_1,
                        MATRIX_3D*  st_MX_2,
                        MATRIX_2D*  sp_MX_2 )
{
   int cmp;
   cmp = MATRIX_3D_Compare( st_MX_1, st_MX_2 );
   if ( cmp != 0 ) return cmp;
   cmp = MATRIX_2D_Compare( sp_MX_1, sp_MX_2 );
   if ( cmp != 0 ) return cmp;
   return 0;
}

/*
 *  FUNCTION:  DP_MATRIX_Diff()
 *  SYNOPSIS: Compare two dynamic programming matrices.
 */
int DP_MATRIX_Diff ( MATRIX_3D*  st_MX_1,
                     MATRIX_2D*  sp_MX_1,
                     MATRIX_3D*  st_MX_2,
                     MATRIX_2D*  sp_MX_2,
                     MATRIX_3D*  st_MX_res,
                     MATRIX_2D*  sp_MX_res )
{
   MATRIX_3D_Diff( st_MX_1, st_MX_2, st_MX_res );
   MATRIX_2D_Diff( sp_MX_1, sp_MX_2, sp_MX_res );
}

/*
 *  FUNCTION:  DP_MATRIX_Save()
 *  SYNOPSIS:  Save dynamic programming matrix to file (by filename).
 */
void DP_MATRIX_Save( const int         Q,
                     const int         T,
                     MATRIX_3D*        st_MX,
                     MATRIX_2D*        sp_MX,
                     const char*       _filename_ )
{
   printf("Saving matrix...\n");
   FILE *fp;
   fp = fopen(_filename_, "w");
   DP_MATRIX_Dump(Q, T, st_MX, sp_MX, fp);
   fclose(fp);
   printf("Saved DP_MATRIX to: '%s'\n", _filename_);
}

/*
 *  FUNCTION:  DP_MATRIX_Dump()
 *  SYNOPSIS:  Output dynamic programming matrix to file.
 */
void DP_MATRIX_Dump( const int         Q,
                     const int         T,
                     MATRIX_3D*        st_MX,
                     MATRIX_2D*        sp_MX,
                     FILE*             fp )
{
   /* padding between numbers */
   int pad     = 9;
   int r_pad   = 4;

   /* PRINT resulting dp matrix */
   fprintf(fp, "##### DP MATRIX ##### \n");
   fprintf(fp, "XDIM\t%d\t%d\n\n", Q, T);

   /* Header */
   fprintf(fp, "##### NORMAL STATES #####\n");
   fprintf(fp, "XMATRIX\n");
   /* Header Indices */
   fprintf(fp, "%*s   ", -3, "#");
   for (int i = 0; i <= T; i++)
   {
      fprintf(fp, "%9d ", i);
   }
   fprintf(fp, "\n");

   /* Row-by-Row Values */
   for (int i = 0; i <= Q; i++)
   {
      fprintf(fp, "%*d M ", -r_pad, i);
      for (int j = 0; j <= T; j++)
      {
         fprintf(fp, "%9.4f ", MMX(i, j) );
      }
      fprintf(fp, "\n");

      fprintf(fp, "%*s I ", -r_pad, "");
      for (int j = 0; j <= T; j++)
      {
         fprintf(fp, "%9.4f ", IMX(i, j) );
      }
      fprintf(fp, "\n");

      fprintf(fp, "%*s D ", -r_pad, "");
      for (int j = 0; j <= T; j++)
      {
         fprintf(fp, "%9.4f ", DMX(i, j) );
      }
      fprintf(fp, "\n\n");
   }
   fprintf(fp, "/\n\n");

   fprintf(fp, "###### SPECIAL STATES #####\n");
   fprintf(fp, "N\t");
   for (int i = 0; i <= Q; i++)
   { fprintf(fp, "%9.4f ", XMX(SP_N, i) ); }
   fprintf(fp, "\n");
   fprintf(fp, "J\t");
   for (int i = 0; i <= Q; i++)
   { fprintf(fp, "%9.4f ", XMX(SP_J, i) ); }
   fprintf(fp, "\n");
   fprintf(fp, "E\t");
   for (int i = 0; i <= Q; i++)
   { fprintf(fp, "%9.4f ", XMX(SP_E, i) ); }
   fprintf(fp, "\n");
   fprintf(fp, "C\t");
   for (int i = 0; i <= Q; i++)
   { fprintf(fp, "%9.4f ", XMX(SP_C, i) ); }
   fprintf(fp, "\n");
   fprintf(fp, "B\t");
   for (int i = 0; i <= Q; i++)
   { fprintf(fp, "%9.4f ", XMX(SP_B, i) ); }
   fprintf(fp, "\n");
}

/*
 *  FUNCTION:  DP_MATRIX_Trace_Save()
 *  SYNOPSIS:  Save dynamic programming matrix with trace to file.
 */
void DP_MATRIX_Trace_Save( const int         Q,
                           const int         T,
                           MATRIX_3D*        st_MX,
                           MATRIX_2D*        sp_MX,
                           ALIGNMENT*        tr,
                           const char*       _filename_ )
{
   FILE *fp;
   fp = fopen(_filename_, "w");
   DP_MATRIX_Trace_Dump( Q, T, st_MX, sp_MX, tr, fp );
   fclose( fp );
   printf("Saved DP_MATRIX with TRACE to: '%s'\n", _filename_);
}


/*
 *  FUNCTION:  DP_MATRIX_Trace_Dump()
 *  SYNOPSIS:  Save dynamic programming matrix with trace to file.
 */
void DP_MATRIX_Trace_Dump( const int         Q,
                           const int         T,
                           MATRIX_3D*        st_MX,
                           MATRIX_2D*        sp_MX,
                           ALIGNMENT*        tr,
                           FILE*             fp )
{
   /* PRINT resulting dp matrix */
   fprintf(fp, "##### DP MATRIX ##### \n");
   fprintf(fp, "X_DIM\t%d\t%d\n\n", Q, T);

   /* Traceback */
   fprintf(fp, "X_TRACE\n");
   /* Traceback */
   for (int i = 0; i < tr->N; i++)
   {
      int st = tr->traces[i].st;
      if ( st == M_ST ) {
         fprintf(fp, "[%d]%s\t%d\t%d\n", i, STATE_NAMES[st], tr->traces[i].i, tr->traces[i].j);
      }
   }
   fprintf(fp, "/\n\n");

   /* Header */
   fprintf(fp, "##### NORMAL STATES #####\n");
   fprintf(fp, "X_MATRIX\n");
   /* Header Indices */
   fprintf(fp, "#\t");
   for (int i = 0; i <= T; i++)
   {
      fprintf(fp, "%d\t", i);
   }
   fprintf(fp, "\n");

   /* Row-by-Row Values */
   for (int i = 0; i <= Q; i++)
   {
      fprintf(fp, "M %d\t", i );
      for (int j = 0; j <= T; j++)
      {
         fprintf(fp, "%.3f\t", MMX(i, j) );
      }
      fprintf(fp, "\n");

      fprintf(fp, "I %d\t", i );
      for (int j = 0; j <= T; j++)
      {
         fprintf(fp, "%.3f\t", IMX(i, j) );
      }
      fprintf(fp, "\n");

      fprintf(fp, "D %d\t", i );
      for (int j = 0; j <= T; j++)
      {
         fprintf(fp, "%.3f\t", DMX(i, j) );
      }
      fprintf(fp, "\n\n");
   }
   fprintf(fp, "/\n\n");

   fprintf(fp, "###### SPECIAL STATES #####\n");
   fprintf(fp, "X_SPECIAL\n");
      /* Header Indices */
   fprintf(fp, "#\t");
   for (int i = 0; i <= Q; i++)
   {
      fprintf(fp, "%d\t", i);
   }
   fprintf(fp, "\n");
   fprintf(fp, "N\t");
   for (int i = 0; i <= Q; i++)
   {
      fprintf(fp, "%.3f\t", XMX(SP_N, i) );
   }
   fprintf(fp, "\n");
   fprintf(fp, "J\t");
   for (int i = 0; i <= Q; i++)
   {
      fprintf(fp, "%.3f\t", XMX(SP_J, i) );
   }
   fprintf(fp, "\n");
   fprintf(fp, "E\t");
   for (int i = 0; i <= Q; i++)
   {
      fprintf(fp, "%.3f\t", XMX(SP_E, i) );
   }
   fprintf(fp, "\n");
   fprintf(fp, "C\t");
   for (int i = 0; i <= Q; i++)
   {
      fprintf(fp, "%.3f\t", XMX(SP_C, i) );
   }
   fprintf(fp, "\n");
   fprintf(fp, "B\t");
   for (int i = 0; i <= Q; i++)
   {
      fprintf(fp, "%.3f\t", XMX(SP_B, i) );
   }
   fprintf(fp, "\n/\n");
}

/*
 *  FUNCTION:  DP_MATRIX_VIZ_Compare()
 *  SYNOPSIS:  Projects two EDGEBOUNDS onto 2D_MATRIX.
 */
void DP_MATRIX_VIZ_Compare(   MATRIX_2D*        cloud_MX,
                              EDGEBOUNDS*       edg_1,
                              EDGEBOUNDS*       edg_2 )
{
   MATRIX_2D_Fill( cloud_MX, 0 );
   MATRIX_2D_Cloud_Fill( cloud_MX, edg_1, 1 );
   MATRIX_2D_Cloud_Fill( cloud_MX, edg_2, 2 );
}

/*
 *  FUNCTION:  DP_MATRIX_VIZ_Trace()
 *  SYNOPSIS:  Adds trace to visulization of dp matrix.
 */
void DP_MATRIX_VIZ_Trace(  MATRIX_2D*        cloud_MX,
                           const ALIGNMENT*  aln )
{
   for ( int i = 0; i < aln->N; i++ ) {
      TRACE* tr = &(aln->traces[i]);
      if ( tr->st == M_ST )
         MX_2D( cloud_MX, tr->i, tr->j ) = -1.0;
   }
}

/*
 *  FUNCTION:  DP_MATRIX_VIZ_Save()
 *  SYNOPSIS:  Saves simple visualization to filename.
 */
void DP_MATRIX_VIZ_Save(   MATRIX_2D*        cloud_MX,
                           char*             filename )
{
   FILE* fp = fopen( filename, "w+" );
   if ( fp == NULL ) {
      
   }
}

/*
 *  FUNCTION:  DP_MATRIX_VIZ_Dump()
 *  SYNOPSIS:  Outputs simple visualization of dp matrix.
 */
void DP_MATRIX_VIZ_Dump(   MATRIX_2D*        cloud_MX,
                           FILE*             fp )
{
   /* padding between points */
   int   pad   = 2;
   /* spacing between labels */
   int   lbl   = 10;
   /* left edge padding for labels */
   int   l_pad = 5;

   /* shaded square

   /* symbol lookup table */
   int   sym_tbl_size   = 6;
   char  syms[6]        = { 'o',  '.', '/',  '\\', 'X', '#' };
   // char  syms[5]        = { 'o',  '.', '1',  '2', '3', '4' };
   // char  syms[5]        = { 'o',  '.', 176,  177, 178 };
   int   vals[6]        = { -1.0, 0.0, 1.0,  2.0, 3.0, 4.0 };

   fprintf(fp, "%*s", -l_pad, "");
   for ( int j = 0; j < cloud_MX->C; j += lbl  ) {
         fprintf(fp, "%*d", -1 * (pad * lbl), j);
   }
   fprintf(fp, "\n");
   

   for ( int i = 0; i < cloud_MX->R; i++ ) {
      if ( i % lbl == 0 )
         fprintf(fp, "%*d", -l_pad, i);
      else
         fprintf(fp, "%*s", -l_pad, "");
      for ( int j = 0; j < cloud_MX->C; j++ ) {
         for ( int k = 0; k < sym_tbl_size; k++ ) {
            if ( MX_2D( cloud_MX, i, j ) == vals[k] )
               fprintf(fp, "%*c", pad, syms[k]);
         }
      }
      fprintf(fp, "\n");
   }
   fprintf(fp, "\n");
}

/*
 *  FUNCTION:  DP_MATRIX_VIZ_Dump()
 *  SYNOPSIS:  Outputs simple visualization of dp matrix.
 *             Symbols are color-coded.
 *             Only works when piped to stdout or stderr.
 */
void DP_MATRIX_VIZ_Color_Dump(   MATRIX_2D*        cloud_MX,
                                 FILE*             fp )
{
   /* if not being piped to stdout or stderr, don't use color */
   if ( fp != stdout && fp != stderr ) {
      DP_MATRIX_VIZ_Dump( cloud_MX, fp );
      return;
   }

   /* padding between points */
   int   pad   = 2;
   /* spacing between labels */
   int   lbl   = 10;
   /* right edge padding for labels */
   int   l_pad = 5;

   /* symbol lookup table */
   int   sym_tbl_size   = 5;
   int   vals[5]        = { -1.0,  0.0,  1.0,  2.0, 3.0 };        /* cloud value */
   char  syms[5]        = { 'o',  '.', '/',  '\\', 'X' };         /* symbol representation of value */
   // char  syms[5]        = { 'o',  '.', 176,  177, 178 };       /* unicode symbols? */
   int   colors[5]      = {    1,    6,    3,    4,   5 };        /* color of symbol */
   bool  bold[5]        = {    1,    0,    1,    1,   1 };        /* whether color is bold or not */

   fprintf(fp, "%*s", -l_pad, "");
   for ( int j = 0; j < cloud_MX->C; j += lbl  ) {
         fprintf(fp, "%*d", -1 * (pad * lbl), j);
   }
   fprintf(fp, "\n");

   for ( int i = 0; i < cloud_MX->R; i++ ) {
      if ( i % lbl == 0 )
         fprintf(fp, "%*d", -l_pad, i);
      else
         fprintf(fp, "%*s", -l_pad, "");
      for ( int j = 0; j < cloud_MX->C; j++ ) {
         for ( int k = 0; k < sym_tbl_size; k++ ) {
            if ( MX_2D( cloud_MX, i, j ) == vals[k] ) {
               TEST_set_color_num( colors[k], true );
               fprintf(fp, "%*c", pad, syms[k]);
               TEST_set_color_num( 0, false );
            }
         }
      }
      fprintf(fp, "\n");
   }
   fprintf(fp, "\n");
}

/*
 *  FUNCTION:  DP_MATRIX_VIZ_Dump()
 *  SYNOPSIS:  Outputs match state of matrix.
 */
void DP_MATRIX_MAT_Dump(   int         Q,
                           int         T,
                           MATRIX_3D*  st_MX,
                           FILE*       fp )
{
   /* Header */
   fprintf(fp, "##### MATRIX - MATCH STATE #####\n");
   /* Header Indices */
   fprintf(fp, "#\t");
   for (int i = 0; i <= T; i++)
   {
      fprintf(fp, "%d\t", i);
   }
   fprintf(fp, "\n");

   /* Row-by-Row Values */
   for (int i = 0; i <= Q; i++)
   {
      fprintf(fp, "M %d\t", i );
      for (int j = 0; j <= T; j++)
      {
         fprintf(fp, "%.2f\t", MMX(i, j) );
      }
      fprintf(fp, "\n");
   }
   fprintf(fp, "//\n\n");
}