/*******************************************************************************
 *  FILE:      error_handler.c
 *  PURPOSE:   Functions for handling error codes
 *
 *  AUTHOR:    Dave Rich
 *  BUG:
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <stdlib.h>
#include <execinfo.h>

#include <unistd.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

/* local imports */
#include "structs.h"
#include "objects.h"

/* header */
#include "utilities.h"

void ERRORCHECK_handler(const int    error_code,
                        const char*  _file_,
                        const int    _line_,
                        const char*  _func_,
                        const char*  err_msg )
{
   /* if no message, output default message according to error code */
   if ( err_msg == NULL ) {
      fprintf(stderr, "ERROR (%d): ", error_code);
      switch ( error_code )
      {
         case ERROR_UNKNOWN: {
            fprintf( stderr, "%s\n", "Unknown error." );
         } break;
         case ERROR_MALLOC: {
            fprintf( stderr, "%s\n", "Malloc error." );
         } break;
         case ERROR_REALLOC: {
            fprintf( stderr, "%s\n", "Realloc error." );
         } break;
         default: {
            fprintf( stderr, "%s\n", "Error occurred with invalid error code." );
         }
      }
   }
   /* report custom message */
   else {
      fprintf(stderr, "ERROR: %s\n", err_msg);
   }
   /* report location of error */
   if ( _file_ != NULL && _func_ != NULL ) {
      fprintf(stderr, "ERROR occurred ==> in FILE: \"%s\", at LINE: %d, in FUNC: \"%s\".\n", _file_, _line_, _func_);
   }
   /* terminate program */
   fprintf(stderr, "Program terminated.\n");
   exit(EXIT_FAILURE);
}

/* open file and check for null file pointer */ 
FILE* ERRORCHECK_fopen( char*       filename,
                        const char* permission,
                        const char*  _file_,
                        const int    _line_,
                        const char*  _func_ )
{
   FILE* fp = fopen( filename, permission );
   if (fp == NULL) {
      ERRORCHECK_handler( ERROR_FILE_IO, _file_, _line_, _func_, NULL );
   }
}

/* malloc and check for null pointer */
void* ERRORCHECK_malloc( const int      size,
                         const char*    _file_,
                         const int      _line_,
                         const char*    _func_ )
{
    void* ptr = malloc( size );
    if ( ptr == NULL ) {
        ERRORCHECK_handler( ERROR_MALLOC, _file_, _line_, _func_, NULL );
    }
    return ptr;
}

/* realloc and check for null pointer */
void* ERRORCHECK_realloc( void*         ptr,
                          const int     size,
                          const char*   _file_,
                          const int     _line_,
                          const char*   _func_ )
{
   ptr = realloc( ptr, size );
   if ( ptr == NULL ) {
      ERRORCHECK_handler( ERROR_REALLOC, _file_, _line_, _func_, NULL );
   }
   return ptr;
}

/* free pointer if not already freed, then set to NULL */
void* ERRORCHECK_free(  void*    ptr )
{
   if ( ptr != NULL ) free( ptr );
   ptr = NULL;
   return ptr;
}


/* realloc and check for null pointer */
void* ERRORCHECK_boundscheck( int            idx,
                              int            max,
                              const int      size,
                              const char*    _file_,
                              const int      _line_,
                              const char*    _func_ )
{
   if (idx > max) {
      ERRORCHECK_handler( ERROR_OUT_OF_BOUNDS, _file_, _line_, _func_, NULL );
   }
}

/* error if called function is currently unsupported */
void ERRORCHECK_unsupported_op( const char*     _file_,
                                const int       _line_,
                                const char*     _func_ )
{

}

/* reports an incorrect value encountered in matrix */
void memcheck_error( int row, int col, float mat, float ins, float del )
{
   printf( "#> ERROR: Memory at position (%d,%d) not cleared. Value = ( %9.4f %9.4f %9.4f )\n", 
      row, col, mat, ins, del );
}

/* print a stacktrace in the event of an error */
void ERRORCHECK_stacktrace()
{
   void*   array[10];
   size_t  size;
   char**  strings;
   size_t  i;

   size = backtrace(array, 10);
   strings = backtrace_symbols(array, size);

   fprintf(stderr, "Obtained %ld stack frames.\n", size);

   for (i = 0; i < size; i++) {
      fprintf(stderr, "%s\n", strings[i]);
   }
   free (strings);
}

/* close program and print stacktrace */
void ERRORCHECK_exit( int exit_flag )
{
   ERRORCHECK_stacktrace();
   exit( exit_flag );
}


