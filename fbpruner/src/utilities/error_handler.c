/*******************************************************************************
 *  FILE:      error_handler.c
 *  PURPOSE:   Functions for handling error codes.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:
 *******************************************************************************/

/* === imports === */
#include <stdio.h>
#include <stdlib.h>
#include <execinfo.h>
#include <unistd.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

/* === local imports === */
#include "../objects/structs.h"
#include "../objects/_objects.h"

/* === package header === */
#include "_utilities.h"

/* === private functions === */
/* NONE */

/* === public functions === */
#include "error_handler.h"

void 
ERRORCHECK_handler(  const ERROR_FLAG      error_code,
                     const char*          _file_,
                     const int            _line_,
                     const char*          _func_,
                     const char*          err_msg )
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
         case ERROR_FILE_IO: {
            fprintf( stderr, "%s\n", "File I/O error." );
         }
         default: {
            fprintf( stderr, "%s\n", "Error occurred with invalid error code." );
         }
      }
   }
   /* report custom message */
   else {
      fprintf(stderr, "# ERROR: %s\n", err_msg);
   }
   /* report location of error */
   if ( _file_ != NULL && _func_ != NULL ) {
      fprintf(stderr, "# ERROR occurred ==> in FILE: \"%s\", at LINE: %d, in FUNC: \"%s\".\n", _file_, _line_, _func_);
   }
   /* terminate program */
   fprintf(stderr, "# Program terminated.\n");
   exit(EXIT_FAILURE);
}


inline
FILE* 
ERRORCHECK_fopen( char*          filename,
                  const char*    permission,
                  const char*    _file_,
                  const int      _line_,
                  const char*    _func_ )
{
   FILE* fp = fopen( filename, permission );
   if (fp == NULL) {
      ERRORCHECK_handler( ERROR_FILE_IO, _file_, _line_, _func_, filename );
   }
}


inline
FILE* 
ERRORCHECK_fclose(   FILE*          fp,
                     const char*    _file_,
                     const int      _line_,
                     const char*    _func_ )
{
   /* close file in not null */
   if (fp != NULL) {
      fclose( fp );
   }
   /* else, attempting to close a NULL file pointer */
   else {
      printf("WARNING: Attempted to close a NULL pointer.");
   }
   
   return NULL;
}


inline
void* 
ERRORCHECK_alloc( void*          ptr,
                  const size_t   size,
                  const char*    _file_,
                  const int      _line_,
                  const char*    _func_ )
{
   /* NOTE: realloc behaves like malloc when {ptr} is NULL. */
   ptr = realloc( ptr, size );
   /* if {ptr} is NULL, then we have a memory error */
   if ( ptr == NULL ) {
      printf("ERROR_SIZE: %ld\n", size);
      ERRORCHECK_handler( ERROR_MALLOC, _file_, _line_, _func_, NULL );
   }
   return ptr;
}


inline
void* 
ERRORCHECK_malloc(   const size_t   size,
                     const char*    _file_,
                     const int      _line_,
                     const char*    _func_ )
{
   void* ptr;
   ptr = malloc( size );
   /* if {ptr} is NULL, then we have a memory error */
   if ( ptr == NULL ) {
      printf("ERROR_SIZE: %ld\n", size);
      ERRORCHECK_handler( ERROR_MALLOC, _file_, _line_, _func_, NULL );
   }
   return ptr;
}


inline
void* 
ERRORCHECK_realloc(  void*          ptr,
                     const size_t   size,
                     const char*    _file_,
                     const int      _line_,
                     const char*    _func_ )
{
   ptr = realloc( ptr, size );
   /* if {ptr} is NULL, then we have a memory error */
   if ( ptr == NULL ) {
      printf("ERROR_SIZE: %ld\n", size);
      ERRORCHECK_handler( ERROR_REALLOC, _file_, _line_, _func_, NULL );
   }
   return ptr;
}


inline
void* 
ERRORCHECK_free(  void*          ptr,
                  const char*    _file_,
                  const int      _line_,
                  const char*    _func_ )
{
   if ( ptr != NULL ) {
      free( ptr );
   }

   return NULL;
}


void* 
ERRORCHECK_boundscheck( int            idx,
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


void 
ERRORCHECK_print_location( FILE*          fp,
                           const char*    _file_,
                           const int      _line_,
                           const char*    _func_ )
{
   fprintf( fp, "{ FILE: %s, LINE: %d, FUNCTION: %s }\n", 
      _file_, _line_, _func_ );
}


void 
ERRORCHECK_unsupported_op(    const char*     _file_,
                              const int       _line_,
                              const char*     _func_ )
{

}


void 
ERRORCHECK_memcheck(    int      row, 
                        int      col, 
                        float    mat, 
                        float    ins, 
                        float    del )
{
   printf( "#> ERROR: Memory at position (%d,%d) not cleared. Value = ( %9.4f %9.4f %9.4f )\n",
           row, col, mat, ins, del );
}


void 
ERRORCHECK_stacktrace()
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


void 
ERRORCHECK_exit( int exit_flag )
{
   ERRORCHECK_stacktrace();
   exit( exit_flag );
}


