/*******************************************************************************
 *  FILE:      error_handler.h
 *  PURPOSE:   Functions for handling error calls
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _ERROR_HANDLER_H
#define _ERROR_HANDLER_H

/* realloc and check for null pointer */
void* ERRORCHECK_realloc( void*         ptr,
                          const int     size,
                          const char*   _file_,
                          const int     _line_,
                          const char*   _func_ );

/* open file and check for null file pointer */ 
FILE* ERRORCHECK_fopen( char*       filename,
                        const char* permission,
                        const char*  _file_,
                        const int    _line_,
                        const char*  _func_ );

/* malloc and check for null pointer */
void* ERRORCHECK_malloc( const int      size,
                         const char*    _file_,
                         const int      _line_,
                         const char*    _func_ );

/* realloc and check for null pointer */
void* ERRORCHECK_realloc( void*          ptr,
                          const int      size,
                          const char*    _file_,
                          const int      _line_,
                          const char*    _func_ );

/* free pointer if not already freed, then set to NULL */
void* ERRORCHECK_free(  void*    ptr );

/* print a stacktrace in the event of an error */
void ERRORCHECK_stacktrace();

/* close program and print stacktrace */
void ERRORCHECK_exit( int exit_flag );

void memcheck_error( int row, int col, float mat, float ins, float del );

#endif /* _ERROR_HANDLER_H */