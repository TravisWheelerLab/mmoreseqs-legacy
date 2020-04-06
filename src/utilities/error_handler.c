/*******************************************************************************
 *  FILE:      error_handler.c
 *  PURPOSE:   Functions for handling error codes
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

/* data structures */
#include "objects/structs.h"
#include "utilities/utility.h"

/* debugging methods */
#include "testing.h"

/* header */
#include "error_handler.h"

void ERROR_Handler( const int    error_code,
                    const char*  file,
                    const int    line,
                    const char*  note )
{
   switch( error_code )
   {
      case ERROR_UNKNOWN:
         fprintf( stderr, "%s\n", "ERROR: Unknown error." );
         break;
      case ERROR_MALLOC:
         fprintf( stderr, "%s\n", "ERROR: Malloc error." );
         break;
      default:
         fprintf( stderr, "%s\n", "ERROR: Error occurred with invalid error code." );
   }
   if (file != NULL)
      fprintf(stderr, "ERROR occurred in FILE: \"%s\", at LINE: %d\n)", file, line);
   fprintf(stderr, "[ CODE = %d ]\n", error_code);
   if (note != NULL)
      fprintf(stderr, "NOTE: %s\n", note);

   fprintf(stderr, "Program terminated.\n");
   exit(EXIT_FAILURE);
}

/* check for null pointer */
void ERROR_NullPtrCheck( const void*   data,
                         int*          error_code,
                         const char*   note )
{

}