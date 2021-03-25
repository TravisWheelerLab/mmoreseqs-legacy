/*******************************************************************************
 *  FILE:      report_util.c
 *  PURPOSE:   Utility functions for generating reports.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       - 
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
#include "../objects/structs.h"
#include "../utilities/_utilities.h"
#include "../objects/_objects.h"
#include "../parsers/_parsers.h"
#include "../algs_linear/_algs_linear.h"
#include "../algs_quad/_algs_quad.h"
#include "../algs_naive/_algs_naive.h"

/* header */
#include "_reporting.h"
#include "report_util.h"

/* === UTILTITY FUNCTIONS === */

/*!  FUNCTION:    REPORT_horizontal_rule()
 *   SYNOPSIS:    Print a horizontal rule. 
 */
inline
STATUS_FLAG 
REPORT_horizontal_rule( FILE*    fp )
{
   fprintf( fp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n" );
}

/*!  FUNCTION:    REPORT_hr_size()
 *   SYNOPSIS:    Print a horizontal rule of custom <L> length.
 */
STATUS_FLAG 
REPORT_hr_size(   FILE*          fp, 
                  const int      L )
{
   char* base_hr = "-";

   fprintf( fp, "#");
   for (int i = 0; i < L; i++) {
      fprintf( fp, "%s", base_hr);
   }
   fprintf( fp, "\n");
}

/*!  FUNCTION:    STRING_hr_size()
 *   SYNOPSIS:    Store a horizontal rule of custom <L> length into <str_buf>.
 *                <str_buf> must be of at least length <L+1>.
 */
STATUS_FLAG 
STRING_hr_size(   char*          str_buf, 
                  const int      L )
{
   for (int i = 0; i < L; i++) {
      str_buf[i] = '-';
   }
   str_buf[L] = '\0';
}

/*!  FUNCTION:    REPORT_header()
 *   SYNOPSIS:    Prints field and horizontal rules for <header> list of length <num_fields>.
 *   NOTE:        Only supports title lengths under 127 chars.
 */
STATUS_FLAG 
REPORT_header(    FILE*          fp, 
                  const char*    headers[],
                  const int      num_fields )
{
   char  hr_buf[256];

   /* headers */
   fprintf( fp, "#" );
   for (int i = 0; i < num_fields; i++) {
      fprintf( fp, "%s\t", headers[i] );
   }
   fprintf( fp, "\n" ); 
      
   /* horizontal rules */
   fprintf( fp, "#" );
   for (int i = 0; i < num_fields; i++) {
      STRING_hr_size( hr_buf, strlen(headers[i]) );
      fprintf( fp, "%s\t", hr_buf );
   }
   fprintf( fp, "\n" ); 

   return STATUS_SUCCESS;
}

/*!  FUNCTION:    REPORT_entry()
 *   SYNOPSIS:    Prints tab-delimited of array <data> with length <L> onto line according to its <type>.
 *                If data is type float or double, <sig_digits> determines number of significant digits.
 */
inline
STATUS_FLAG 
REPORT_entry(     FILE*             fp, 
                  const GEN    data[],
                  const int         num_fields,
                  const int         sig_digits )
{
   /* data array */
   for (int i = 0; i < num_fields; i++) {
      REPORT_data(fp, data[i], sig_digits);
   }
   fprintf( fp, "\n" ); 

   return STATUS_SUCCESS;
}

/*!  FUNCTION:    REPORT_horizontal()
 *   SYNOPSIS:    Prints array of <data>, one entry per line, with <header> labels.
 */
STATUS_FLAG 
REPORT_entry_multiline(    FILE*             fp, 
                           const char*       headers[],
                           const GEN    data[],
                           const int         num_fields,
                           const int         sig_digits )
{
   char     str_buf[128];
   size_t   buf_size    = 128;
   int      pad         = 32;

   /* data array */
   for (int i = 0; i < num_fields; i++) {
      GEN_ToString( data[i], str_buf, buf_size, sig_digits );
      fprintf( fp, "%*s:\t%s\n", pad, headers[i], str_buf );
   }
   fprintf( fp, "\n"); 

   return STATUS_SUCCESS;
}

/*!  FUNCTION:    REPORT_data()
 *   SYNOPSIS:    Prints one entry <data> to an open file pointer <fp>.
 */
inline
STATUS_FLAG 
REPORT_data(   FILE*             fp, 
               const GEN    data,
               const int         sig_digits )
{
   char     str_buf[128];
   size_t   buf_size = 128;

   GEN_ToString( data, str_buf, buf_size, sig_digits );
   fprintf( fp, "%s", str_buf);

   return STATUS_SUCCESS;
}