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

#endif /* _ERROR_HANDLER_H */