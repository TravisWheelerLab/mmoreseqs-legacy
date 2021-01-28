/*******************************************************************************
 *  FILE:      report_util.h
 *  PURPOSE:   Reporting Subroutines for generating output.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _REPORT_UTIL_H
#define _REPORT_UTIL_H

/* === macros === */

/* === public functions === */

/*!  FUNCTION:    REPORT_horizontal_rule()
 *   SYNOPSIS:    Print a horizontal rule. 
 */
STATUS_FLAG 
REPORT_horizontal_rule( FILE*    fp );

/*!  FUNCTION:    REPORT_hr_size()
 *   SYNOPSIS:    Print a horizontal rule of specified length.
 */
STATUS_FLAG 
REPORT_hr_size(   FILE*          fp, 
                  const int      L );

/*!  FUNCTION:    STR_hr_size()
 *   SYNOPSIS:    Create a horizonral rule of specified length.
 */
STATUS_FLAG
STRING_hr_size(   char*          str_buf,
                  const int      L );

/*!  FUNCTION:    REPORT_header()
 *   SYNOPSIS:    Prints field headers with 
 */
STATUS_FLAG 
REPORT_header(    FILE*          fp, 
                  const char*    headers[],
                  const int      num_fields );

/*!  FUNCTION:    REPORT_entry()
 *   SYNOPSIS:    Prints tab-delimited of array <data> with length <L> onto line according to its <type>.
 *                If data is type float or double, <sig_digits> determines number of significant digits.
 */
STATUS_FLAG 
REPORT_entry(     FILE*             fp, 
                  const GEN_DATA    data[],
                  const int         num_fields,
                  const int         sig_digits );

/*!  FUNCTION:    REPORT_horizontal()
 *   SYNOPSIS:    Prints array of <data>, one entry per line, with <header> labels.
 */
STATUS_FLAG 
REPORT_entry_multiline(    FILE*             fp, 
                           const char*       headers[],
                           const GEN_DATA    data[],
                           const int         num_fields,
                           const int         sig_digits );

/*!  FUNCTION:    REPORT_data()
 *   SYNOPSIS:    Prints one entry <data> to an open file pointer <fp>.
 */
STATUS_FLAG 
REPORT_data(   FILE*             fp, 
               const GEN_DATA    data,
               const int         sig_digits );


#endif /* _REPORT_UTIL_H */