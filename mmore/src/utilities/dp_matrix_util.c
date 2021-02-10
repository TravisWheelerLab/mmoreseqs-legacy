/*******************************************************************************
 *  FILE:      dp_matrix.h
 *  SYNOPSIS:  Dynamic Programming Matrix functions.
 *
 *  AUTHOR:    Dave Rich
 *  TODO:      This should all be moved into the DP_MATRIX data type when it is finished.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/* local imports */
#include "../objects/structs.h"
#include "../objects/_objects.h"

/* header */
#include "_utilities.h"
#include "dp_matrix_util.h"

/* TODO:*/
/*! FUNCTION:  DP_MATRIX_GetBounds()
 *  SYNOPSIS:  Get the edgebounds of matrix at given antidiagonal (closed form).
 *             Stores new bounds in BOUND <bnd>.
 */
inline
void 
DP_MATRIX_GetBounds(   const int   Q,
                        const int   T,
                        int         d_0,
                        int         d_cnt,
                        BOUND*      bnd )
{
   BOUND b;
}

/*! FUNCTION:  DP_MATRIX_Copy()
 *  SYNOPSIS:  Copy dynamic programming matrix into destination.
 */
void 
DP_MATRIX_Copy(   const int   Q,
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

/*! FUNCTION:  DP_MATRIX_Fill()
 *  SYNOPSIS:  Fill entire dynamic programming matrix with value
 */
void 
DP_MATRIX_Fill(   const int   Q, 
                  const int   T,
                  MATRIX_3D*  st_MX,
                  MATRIX_2D*  sp_MX,
                  const float val )
{
   MATRIX_3D_Fill( st_MX, val );
   MATRIX_2D_Fill( sp_MX, val );
}

/*! FUNCTION:  DP_MATRIX_Clean()
 *  SYNOPSIS:  Fill entire dynamic programming matrix with clean value (-INF).
 */
void 
DP_MATRIX_Clean(  const int   Q, 
                  const int   T,
                  MATRIX_3D*  st_MX,
                  MATRIX_2D*  sp_MX )
{
   if (st_MX != NULL) MATRIX_3D_Clean( st_MX );
   if (sp_MX != NULL) MATRIX_2D_Clean( sp_MX );
}

/*! FUNCTION:  DP_MATRIX_Clean_Verify()
 *  SYNOPSIS:  Check whether there are clean.  If clean, returns true.
 */
bool 
DP_MATRIX_Clean_Verify(    const int   Q, 
                           const int   T,
                           MATRIX_3D*  st_MX,
                           MATRIX_2D*  sp_MX )
{
   return st_MX->clean && sp_MX->clean;
}

/*! FUNCTION:  DP_MATRIX_Save()
 *  SYNOPSIS:  Compare two dynamic programming matrices. 
 *    RETURN:  Returns 0 if equal; otherwise returns count of differing cells.
 */
int 
DP_MATRIX_Compare(   MATRIX_3D*  st_MX_1,
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

/*! FUNCTION:  DP_MATRIX_Diff()
 *  SYNOPSIS: Compare two dynamic programming matrices.
 */
int 
DP_MATRIX_Diff(   MATRIX_3D*  st_MX_1,
                  MATRIX_2D*  sp_MX_1,
                  MATRIX_3D*  st_MX_2,
                  MATRIX_2D*  sp_MX_2,
                  MATRIX_3D*  st_MX_res,
                  MATRIX_2D*  sp_MX_res )
{
   MATRIX_3D_Diff( st_MX_1, st_MX_2, st_MX_res );
   MATRIX_2D_Diff( sp_MX_1, sp_MX_2, sp_MX_res );
}

/*! FUNCTION:  DP_MATRIX_Add()
 *  SYNOPSIS: Add two dynamic programming matrices.
 */
int 
DP_MATRIX_Add(    MATRIX_3D*  st_MX_1,
                  MATRIX_2D*  sp_MX_1,
                  MATRIX_3D*  st_MX_2,
                  MATRIX_2D*  sp_MX_2,
                  MATRIX_3D*  st_MX_res,
                  MATRIX_2D*  sp_MX_res )
{
   MATRIX_3D_Add( st_MX_1, st_MX_2, st_MX_res );
   MATRIX_2D_Add( sp_MX_1, sp_MX_2, sp_MX_res );
}

/*! FUNCTION:  DP_MATRIX_Save()
 *  SYNOPSIS:  Save dynamic programming matrix to file (by filename).
 */
void 
DP_MATRIX_Save(   const int         Q,
                  const int         T,
                  MATRIX_3D*        st_MX,
                  MATRIX_2D*        sp_MX,
                  const char*       filename )
{
   printf("Saving matrix...\n");
   FILE *fp;
   fp = fopen(filename, "w");
   DP_MATRIX_Dump(Q, T, st_MX, sp_MX, fp);
   fclose(fp);
   printf("Saved DP_MATRIX to: '%s'\n", filename);
}

/*! FUNCTION:  DP_MATRIX_Dump()
 *  SYNOPSIS:  Output dynamic programming matrix to file.
 */
void 
DP_MATRIX_Dump(   const int         Q,
                  const int         T,
                  MATRIX_3D*        st_MX,
                  MATRIX_2D*        sp_MX,
                  FILE*             fp )
{
   int pad, dec;
   pad = 9;
   dec = 3;

   /* PRINT resulting dp matrix */
   fprintf(fp, "#####_MMORE_DP_MATRIX_##### \n");
   fprintf(fp, "XDIM\t%d\t%d\n\n", Q, T);

   /* Header */
   fprintf(fp, "%*s ", -14, "# #");
   for (int i = 0; i <= T; i++) {
      fprintf(fp, "%*d ", -9, i);
   }
   fprintf(fp, "\n");

   /* Row-by-Row */
   for (int i = 0; i <= Q; i++)
   {
      fprintf(fp, "M %*d ", -4, i );
      for (int j = 0; j <= T; j++) {
         fprintf(fp, "%*.*f ", pad, dec, MMX(i, j) );
      }
      fprintf(fp, "\n");

      fprintf(fp, "I %*d ", -4, i );
      for (int j = 0; j <= T; j++) {
         fprintf(fp, "%*.*f ", pad, dec, IMX(i, j) );
      }
      fprintf(fp, "\n");

      fprintf(fp, "D %*d ", -4, i );
      for (int j = 0; j <= T; j++) {
         fprintf(fp, "%*.*f ", pad, dec, DMX(i, j) );
      }
      fprintf(fp, "\n\n");
   }

   fprintf(fp, "###### SPECIAL STATES #####\n");
   fprintf(fp, "%*s ", -14, "#");
   for (int i = 0; i <= Q; i++) {
      fprintf(fp, "%*d ", -9, i);
   }
   fprintf(fp, "\n");

   fprintf(fp, "%*s ", -6, "N");
   for (int i = 0; i <= Q; i++) { 
      fprintf(fp, "%*.*f ", pad, dec, XMX(SP_N, i) ); 
   }
   fprintf(fp, "\n");
   
   fprintf(fp, "%*s ", -6, "J");
   for (int i = 0; i <= Q; i++) { 
      fprintf(fp, "%*.*f ", pad, dec, XMX(SP_J, i) ); 
   }
   fprintf(fp, "\n");

   fprintf(fp, "%*s ", -6, "E");
   for (int i = 0; i <= Q; i++) { 
      fprintf(fp, "%*.*f ", pad, dec, XMX(SP_E, i) ); 
   }
   fprintf(fp, "\n");

   fprintf(fp, "%*s ", -6, "C");
   for (int i = 0; i <= Q; i++) { 
      fprintf(fp, "%*.*f ", pad, dec, XMX(SP_C, i) ); 
   }
   fprintf(fp, "\n");

   fprintf(fp, "%*s ", -6, "B");
   for (int i = 0; i <= Q; i++) { 
      fprintf(fp, "%*.*f ", pad, dec, XMX(SP_B, i) ); 
   }
   fprintf(fp, "\n\n");
}


/*! FUNCTION:  DP_MATRIX_Log()
 *  SYNOPSIS:  Convert all dp matrix cells to log space
 */
void 
DP_MATRIX_Log(    const int         Q,
                  const int         T,
                  MATRIX_3D*        st_MX,
                  MATRIX_2D*        sp_MX )
{
   if (sp_MX != NULL) MATRIX_2D_Log( sp_MX );
   if (st_MX != NULL) MATRIX_3D_Log( st_MX );
}


/*! FUNCTION:  DP_MATRIX_Exp()
 *  SYNOPSIS:  Convert all dp matrix cells to normal space
 */
void 
DP_MATRIX_Exp(    const int         Q,
                  const int         T,
                  MATRIX_3D*        st_MX,
                  MATRIX_2D*        sp_MX )
{
   if (sp_MX != NULL) MATRIX_2D_Log( sp_MX );
   if (st_MX != NULL) MATRIX_3D_Log( st_MX );
}


/*! FUNCTION:  DP_MATRIX_Log_Dump()
 *  SYNOPSIS:  Output dynamic programming matrix to file.
 */
void 
DP_MATRIX_Log_Dump(  const int         Q,
                     const int         T,
                     MATRIX_3D*        st_MX,
                     MATRIX_2D*        sp_MX,
                     FILE*             fp )
{
   /* PRINT resulting dp matrix */
   int pad, dec;
   pad = 9;
   dec = 3;

   fprintf(fp, "##### DP MATRIX ##### \n");
   fprintf(fp, "XDIM\t%d\t%d\n\n", Q, T);

   /* Header */
   fprintf(fp, "%*s ", -14, "# #");
   for (int i = 0; i <= T; i++) {
      fprintf(fp, "%*d ", -9, i);
   }
   fprintf(fp, "\n");

   /* Row-by-Row */
   for (int i = 0; i <= Q; i++)
   {
      fprintf(fp, "M %*d ", -4, i );
      for (int j = 0; j <= T; j++) {
         fprintf(fp, "%*.*f ", pad, dec, log( MMX(i, j) ) );
      }
      fprintf(fp, "\n");

      fprintf(fp, "I %*d ", -4, i );
      for (int j = 0; j <= T; j++) {
         fprintf(fp, "%*.*f ", pad, dec, log( IMX(i, j) ) );
      }
      fprintf(fp, "\n");

      fprintf(fp, "D %*d ", -4, i );
      for (int j = 0; j <= T; j++) {
         fprintf(fp, "%*.*f ", pad, dec, log( DMX(i, j) ) );
      }
      fprintf(fp, "\n\n");
   }

   fprintf(fp, "###### SPECIAL STATES #####\n");
   fprintf(fp, "%*s ", -14, "#");
   for (int i = 0; i <= Q; i++) {
      fprintf(fp, "%*d ", -9, i);
   }
   fprintf(fp, "\n");

   fprintf(fp, "%*s ", -6, "N");
   for (int i = 0; i <= Q; i++) { 
      fprintf(fp, "%*.*f ", pad, dec, log( XMX(SP_N, i) ) ); 
   }
   fprintf(fp, "\n");
   
   fprintf(fp, "%*s ", -6, "J");
   for (int i = 0; i <= Q; i++) { 
      fprintf(fp, "%*.*f ", pad, dec, log( XMX(SP_J, i) ) ); 
   }
   fprintf(fp, "\n");

   fprintf(fp, "%*s ", -6, "E");
   for (int i = 0; i <= Q; i++) { 
      fprintf(fp, "%*.*f ", pad, dec, log( XMX(SP_E, i) ) ); 
   }
   fprintf(fp, "\n");

   fprintf(fp, "%*s ", -6, "C");
   for (int i = 0; i <= Q; i++) { 
      fprintf(fp, "%*.*f ", pad, dec, log( XMX(SP_C, i) ) ); 
   }
   fprintf(fp, "\n");

   fprintf(fp, "%*s ", -6, "B");
   for (int i = 0; i <= Q; i++) { 
      fprintf(fp, "%*.*f ", pad, dec, log( XMX(SP_B, i) ) ); 
   }
   fprintf(fp, "\n\n");
}

/** FUNCTION:  DP_MATRIX_Norm_Dump()
 *  SYNOPSIS:  Output dynamic programming normal state matrix {M,I,D} to file.
 */
void 
DP_MATRIX_Norm_Dump(    const int         Q,
                        const int         T,
                        MATRIX_3D*        st_MX,
                        MATRIX_2D*        sp_MX,
                        FILE*             fp )
{
   /* PRINT resulting dp matrix */
   fprintf(fp, "##### DP MATRIX ##### \n");
   fprintf(fp, "XDIM\t%d\t%d\n\n", Q, T);

   /* Header */
   fprintf(fp, "%*s ", -9, "#");
   for (int i = 0; i <= T; i++) {
      fprintf(fp, "%*d ", -7, i);
   }
   fprintf(fp, "\n");

   /* Row-by-Row */
   for (int i = 0; i <= Q; i++)
   {
      fprintf(fp, "M %*d ", -4, i );
      for (int j = 0; j <= T; j++) {
         fprintf(fp, "%7.3f ", MMX(i, j) );
      }
      fprintf(fp, "\n");

      fprintf(fp, "I %*d ", -4, i );
      for (int j = 0; j <= T; j++) {
         fprintf(fp, "%7.3f ", IMX(i, j) );
      }
      fprintf(fp, "\n");

      fprintf(fp, "D %*d ", -4, i );
      for (int j = 0; j <= T; j++) {
         fprintf(fp, "%7.3f ", DMX(i, j) );
      }
      fprintf(fp, "\n\n");
   }
}

/*! FUNCTION:  DP_MATRIX_Dump()
 *  SYNOPSIS:  Output dynamic programming special state matrix to file.
 */
void 
DP_MATRIX_Spec_Dump(    const int         Q,
                        const int         T,
                        MATRIX_3D*        st_MX,
                        MATRIX_2D*        sp_MX,
                        FILE*             fp )
{
   fprintf(fp, "###### SPECIAL STATES #####\n");

   fprintf(fp, "%*s ", -9, "#");
   for (int i = 0; i <= Q; i++) {
      fprintf(fp, "%*d ", -7, i);
   }
   fprintf(fp, "\n");

   fprintf(fp, "%*s ", -6, "N");
   for (int i = 0; i <= Q; i++) { 
      fprintf(fp, "%7.3f ", XMX(SP_N, i) ); 
   }
   fprintf(fp, "\n");
   
   fprintf(fp, "%*s ", -6, "J");
   for (int i = 0; i <= Q; i++) { 
      fprintf(fp, "%7.3f ", XMX(SP_J, i) ); 
   }
   fprintf(fp, "\n");

   fprintf(fp, "%*s ", -6, "E");
   for (int i = 0; i <= Q; i++) { 
      fprintf(fp, "%7.3f ", XMX(SP_E, i) ); 
   }
   fprintf(fp, "\n");

   fprintf(fp, "%*s ", -6, "C");
   for (int i = 0; i <= Q; i++) { 
      fprintf(fp, "%7.3f ", XMX(SP_C, i) ); 
   }
   fprintf(fp, "\n");

   fprintf(fp, "%*s ", -6, "B");
   for (int i = 0; i <= Q; i++) { 
      fprintf(fp, "%7.3f ", XMX(SP_B, i) ); 
   }
   fprintf(fp, "\n");
}

/*! FUNCTION:  DP_MATRIX_Dump()
 *  SYNOPSIS:  Output dynamic programming matrix to file.
 */
void 
DP_MATRIX_SpecExp_Dump(    const int         Q,
                           const int         T,
                           MATRIX_3D*        st_MX,
                           MATRIX_2D*        sp_MX,
                           FILE*             fp )
{
   fprintf(fp, "###### SPECIAL STATES #####\n");

   fprintf(fp, "%*s ", -9, "#");
   for (int i = 0; i <= Q; i++) {
      fprintf(fp, "%*d ", -7, i);
   }
   fprintf(fp, "\n");

   fprintf(fp, "%*s ", -6, "N");
   for (int i = 0; i <= Q; i++) { 
      fprintf(fp, "%7.3f ", exp(XMX(SP_N, i)) ); 
   }
   fprintf(fp, "\n");
   
   fprintf(fp, "%*s ", -6, "J");
   for (int i = 0; i <= Q; i++) { 
      fprintf(fp, "%7.3f ", exp(XMX(SP_J, i)) ); 
   }
   fprintf(fp, "\n");

   fprintf(fp, "%*s ", -6, "E");
   for (int i = 0; i <= Q; i++) { 
      fprintf(fp, "%7.3f ", exp(XMX(SP_E, i)) ); 
   }
   fprintf(fp, "\n");

   fprintf(fp, "%*s ", -6, "C");
   for (int i = 0; i <= Q; i++) { 
      fprintf(fp, "%7.3f ", exp(XMX(SP_C, i)) ); 
   }
   fprintf(fp, "\n");

   fprintf(fp, "%*s ", -6, "B");
   for (int i = 0; i <= Q; i++) { 
      fprintf(fp, "%7.3f ", exp(XMX(SP_B, i)) ); 
   }
   fprintf(fp, "\n");
}

/*! FUNCTION:  DP_MATRIX_Trace_Save()
 *  SYNOPSIS:  Save dynamic programming matrix with trace to file.
 */
void 
DP_MATRIX_Trace_Save(   const int         Q,
                        const int         T,
                        MATRIX_3D*        st_MX,
                        MATRIX_2D*        sp_MX,
                        ALIGNMENT*        tr,
                        const char*       filename )
{
   FILE *fp;
   fp = fopen(filename, "w");
   DP_MATRIX_Trace_Dump( Q, T, st_MX, sp_MX, tr, fp );
   fclose( fp );
   printf("Saved DP_MATRIX with TRACE to: '%s'\n", filename);
}


/*! FUNCTION:  DP_MATRIX_Trace_Dump()
 *  SYNOPSIS:  Save dynamic programming matrix with trace to file.
 */
void 
DP_MATRIX_Trace_Dump(   const int         Q,
                        const int         T,
                        MATRIX_3D*        st_MX,
                        MATRIX_2D*        sp_MX,
                        ALIGNMENT*        tr,
                        FILE*             fp )
{
   int N;
   TRACE* trace;

   N     = tr->traces->N;
   trace = tr->traces->data;

   /* PRINT resulting dp matrix */
   fprintf(fp, "##### DP MATRIX ##### \n");
   fprintf(fp, "X_DIM\t%d\t%d\n\n", Q, T);

   /* Traceback */
   fprintf(fp, "X_TRACE\n");
   /* Traceback */
   for (int i = 0; i < N; i++)
   {
      int st = trace[i].st;
      if ( st == M_ST ) {
         fprintf(fp, "[%d]%s\t%d\t%d\n", i, STATE_NAMES[st], trace[i].q_0, trace[i].t_0);
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

/*! FUNCTION:  DP_MATRIX_VIZ_Compare()
 *  SYNOPSIS:  Projects two EDGEBOUNDS onto 2D_MATRIX.
 */
void 
DP_MATRIX_VIZ_Compare(     MATRIX_2D*        cloud_MX,
                           EDGEBOUNDS*       edg_1,
                           EDGEBOUNDS*       edg_2 )
{
   MATRIX_2D_Fill( cloud_MX, 0.0 );
   MATRIX_2D_Cloud_Fill( cloud_MX, edg_1, 1.0 );
   MATRIX_2D_Cloud_Fill( cloud_MX, edg_2, 2.0 );
}

/*! FUNCTION:  DP_MATRIX_VIZ_Trace()
 *  SYNOPSIS:  Adds trace to visulization of dp matrix.
 */
void 
DP_MATRIX_VIZ_Trace(    MATRIX_2D*        cloud_MX,
                        const ALIGNMENT*  aln )
{
   for ( int i = 0; i < aln->traces->N; i++ ) {
      TRACE* tr = &(aln->traces->data[i]);
      if ( tr->st == M_ST )
         MX_2D( cloud_MX, tr->q_0, tr->t_0 ) = -1.0;
   }
}

/*! FUNCTION:  DP_MATRIX_VIZ_Save()
 *  SYNOPSIS:  Saves simple visualization to filename.
 */
void 
DP_MATRIX_VIZ_Save(     MATRIX_2D*        cloud_MX,
                        char*             filename )
{
   FILE* fp = fopen( filename, "w+" );
   if ( fp == NULL ) {
      
   }
}

/*! FUNCTION:  DP_MATRIX_VIZ_Dump()
 *  SYNOPSIS:  Outputs simple visualization of dp matrix.
 */
void 
DP_MATRIX_VIZ_Dump(     MATRIX_2D*        cloud_MX,
                        FILE*             fp )
{
   /* padding between points */
   int   pad   = 2;
   /* spacing between labels */
   int   lbl   = 10;
   /* left edge padding for labels */
   int   l_pad = 5;
   /* checks if value is covered by symbols */
   bool  key_found;
   float val;
   char  unknown_sym = '?';

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
         key_found = false;
         val = MX_2D( cloud_MX, i, j );
         for ( int k = 0; k < sym_tbl_size; k++ ) {
            if ( val == vals[k] ) {
               fprintf(fp, "%*c", pad, syms[k]);
               key_found = true;
            }
         }
         if ( key_found == false ) {
            fprintf(fp, "%*c", pad, unknown_sym );
         }
      }
      fprintf(fp, "\n");
   }
   fprintf(fp, "\n");
}

/*! FUNCTION:  DP_MATRIX_VIZ_Dump()
 *  SYNOPSIS:  Outputs simple visualization of dp matrix.
 *             Symbols are color-coded.
 *             Only works when piped to stdout or stderr.
 */
void 
DP_MATRIX_VIZ_Color_Dump(     MATRIX_2D*        cloud_MX,
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

/*! FUNCTION:  DP_MATRIX_MAT_Dump()
 *  SYNOPSIS:  Outputs match states of matrix.
 */
void 
DP_MATRIX_MAT_Dump(     int         Q,
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