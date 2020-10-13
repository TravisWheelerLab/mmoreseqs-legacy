/*******************************************************************************
 *  FILE:      error_handler.c
 *  PURPOSE:   Functions for handling error codes.
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
#include "../objects/structs.h"
#include "../objects/objects.h"

/* header */
#include "utilities.h"


/*
 *  FUNCTION:  ERRORCHECK_handler()
 *  SYNOPSIS:  Passed {error_code} is handled.
 */
void ERRORCHECK_handler(   const int      error_code,
                           const char*    _file_,
                           const int      _line_,
                           const char*    _func_,
                           const char*    err_msg )
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
            fprintf( stderr, "%s %s\n", "File I/O error:", err_msg );
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

/*
 *  FUNCTION:  ERRORCHECK_fopen()
 *  SYNOPSIS:  Opens file and returns file pointer.  
 *             If it returns NULL pointer, then error is thrown.
 *             Handles the error messaging, with program location, and closes program.
 */
inline
FILE* ERRORCHECK_fopen( char*          filename,
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

/*
 *  FUNCTION:  ERRORCHECK_fclose()
 *  SYNOPSIS:  Closes file and returns NULL pointer.
 */
inline
void* ERRORCHECK_fclose(   FILE*          fp,
                           const char*    _file_,
                           const int      _line_,
                           const char*    _func_ )
{
   if (fp != NULL) {
      fclose( fp );
   }
   else {
      /* attempting to close a NULL file pointer */
   }
   
   return NULL;
}

/*
 *  FUNCTION:  ERRORCHECK_alloc()
 *  SYNOPSIS:  Allocates (or reallocates) memory and returns a pointer.
 *             If memory error causes a NULL pointer to be returned, then memory error is thrown.
 *             Handles the error messaging, with program location, and closes program.
 */
inline
void* ERRORCHECK_alloc( void*          ptr,
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

/*
 *  FUNCTION:  ERRORCHECK_malloc()
 *  SYNOPSIS:  Allocates (or reallocates) memory and returns a pointer.
 *             If memory error causes a NULL pointer to be returned, then memory error is thrown.
 *             Handles the error messaging, with program location, and closes program.
 */
inline
void* ERRORCHECK_malloc(   const size_t   size,
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

/*
 *  FUNCTION:  ERRORCHECK_realloc()
 *  SYNOPSIS:  Allocates (or reallocates) memory and returns a pointer.
 *             If memory error causes a NULL pointer to be returned, then memory error is thrown.
 *             Handles the error messaging, with program location, and closes program.
 */
inline
void* ERRORCHECK_realloc(  void*          ptr,
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

/*
 *  FUNCTION:  ERRORCHECK_free()
 *  SYNOPSIS:  Allocates (or reallocates) memory and returns a pointer.
 *             If memory error causes a NULL pointer to be returned, then memory error is thrown.
 *             Handles the error messaging, with program location, and closes program.
 */
inline
void* ERRORCHECK_free(  void*          ptr,
                        const char*    _file_,
                        const int      _line_,
                        const char*    _func_ )
{
   if ( ptr != NULL ) {
      free( ptr );
   }

   return NULL;
}

/*
 *  FUNCTION:  ERRORCHECK_boundscheck()
 *  SYNOPSIS:  
 */
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

/*
 *  FUNCTION:  ERRORCHECK_print_location()
 *  SYNOPSIS:  Prints program location given by ERRORCHECK function.
 */
void ERRORCHECK_print_location(  FILE*          fp,
                                 const char*    _file_,
                                 const int      _line_,
                                 const char*    _func_ )
{
   fprintf( fp, "{ FILE: %s, LINE: %d, FUNCTION: %s }\n", 
      _file_, _line_, _func_ );
}

/*
 *  FUNCTION:  ERRORCHECK_unsupported_op()
 *  SYNOPSIS:  Error called when function is currently unsupported.
 */
void ERRORCHECK_unsupported_op( const char*     _file_,
                                const int       _line_,
                                const char*     _func_ )
{

}

/*
 *  FUNCTION:  ERRORCHECK_memcheck()
 *  SYNOPSIS:  Reports an incorrect value encountered in matrix.
 */
void ERRORCHECK_memcheck(  int      row, 
                           int      col, 
                           float    mat, 
                           float    ins, 
                           float    del )
{
   printf( "#> ERROR: Memory at position (%d,%d) not cleared. Value = ( %9.4f %9.4f %9.4f )\n",
           row, col, mat, ins, del );
}

/*
 *  FUNCTION:  ERRORCHECK_stacktrace()
 *  SYNOPSIS:  Print stacktrace.
 */
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

/*
 *  FUNCTION:  ERRORCHECK_exit()
 *  SYNOPSIS:  Exit program and print stacktrace.
 */
void ERRORCHECK_exit( int exit_flag )
{
   ERRORCHECK_stacktrace();
   exit( exit_flag );
}


