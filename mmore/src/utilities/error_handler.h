/*******************************************************************************
 *  FILE:      error_handler.h
 *  PURPOSE:   Functions for handling errors.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _ERROR_HANDLER_H
#define _ERROR_HANDLER_H

/* === imports === */
/* NONE */

/* === macros === */
/* NONE */

/* === public functions === */

/*! FUNCTION:  ERRORCHECK_handler()
 *  SYNOPSIS:  Passed <error_code> is handled.
 */
void 
ERRORCHECK_handler(  const ERROR_FLAG     error_code,
                     const char*          _file_,
                     const int            _line_,
                     const char*          _func_,
                     const char*          err_msg );

/*! FUNCTION:  ERRORCHECK_fopen()
 *  SYNOPSIS:  Opens file.  If it returns NULL pointer, then error is thrown.
 *             Handles the error messaging, with program location, and closes program.
 */
FILE* 
ERRORCHECK_fopen(    const char*    filename,
                     const char*    permission,
                     const char*    _file_,
                     const int      _line_,
                     const char*    _func_ );

/*! FUNCTION:  ERRORCHECK_fclose()
 *  SYNOPSIS:  Closes file and returns NULL pointer.
 */
FILE* 
ERRORCHECK_fclose(   FILE*          fp,
                     const char*    _file_,
                     const int      _line_,
                     const char*    _func_ );

/*! FUNCTION:  ERRORCHECK_alloc()
 *  SYNOPSIS:  Allocates (or reallocates) memory and returns a pointer.
 *             If memory error causes a NULL pointer to be returned, then memory error is thrown.
 *             Handles the error messaging, with program location, and closes program.
 */
void* 
ERRORCHECK_alloc( void*          ptr,
                  const size_t   size,
                  const char*    _file_,
                  const int      _line_,
                  const char*    _func_ );

/*! FUNCTION:  ERRORCHECK_malloc()
 *  SYNOPSIS:  Allocates (or reallocates) memory and returns a pointer.
 *             If memory error causes a NULL pointer to be returned, then memory error is thrown.
 *             Handles the error messaging, with program location, and closes program.
 */
void* 
ERRORCHECK_malloc(   const size_t   size,
                     const char*    _file_,
                     const int      _line_,
                     const char*    _func_ );

/*! FUNCTION:  ERRORCHECK_malloc()
 *  SYNOPSIS:  Allocates (or reallocates) memory and returns a pointer.
 *             If memory error causes a NULL pointer to be returned, then memory error is thrown.
 *             Handles the error messaging, with program location, and closes program.
 */
void* 
ERRORCHECK_realloc(  void*          ptr,
                     const size_t   size,
                     const char*    _file_,
                     const int      _line_,
                     const char*    _func_ );

/*! FUNCTION:  ERRORCHECK_alloc()
 *  SYNOPSIS:  Allocates (or reallocates) memory and returns a pointer.
 *             If memory error causes a NULL pointer to be returned, then memory error is thrown.
 *             Handles the error messaging, with program location, and closes program.
 */
void* 
ERRORCHECK_free(  void*          ptr,
                  const char*    _file_,
                  const int      _line_,
                  const char*    _func_ );

/*! FUNCTION:  ERRORCHECK_boundscheck()
 *  SYNOPSIS:  
 */
void* 
ERRORCHECK_boundscheck(    int            idx,
                           int            max,
                           const int      size,
                           const char*    _file_,
                           const int      _line_,
                           const char*    _func_ );

/*! FUNCTION:  ERRORCHECK_print_location()
 *  SYNOPSIS:  Prints program location given by ERRORCHECK function.
 */
void 
ERRORCHECK_print_location(    FILE*          fp,
                              const char*    _file_,
                              const int      _line_,
                              const char*    _func_ );

/*! FUNCTION:  ERRORCHECK_unsupported_op()
 *  SYNOPSIS:  Error called when function is currently unsupported.
 */
void ERRORCHECK_unsupported_op( const char*     _file_,
                                const int       _line_,
                                const char*     _func_ );

/*! FUNCTION:  ERRORCHECK_memcheck()
 *  SYNOPSIS:  Reports an incorrect value encountered in matrix.
 */
void ERRORCHECK_memcheck(  int      row, 
                           int      col, 
                           float    mat, 
                           float    ins, 
                           float    del );

/*! FUNCTION:  ERRORCHECK_stacktrace()
 *  SYNOPSIS:  Print stacktrace.
 */
void ERRORCHECK_stacktrace();

/*! FUNCTION:  ERRORCHECK_unsupported_op()
 *  SYNOPSIS:  Error called when function is currently unsupported.
 */
void ERRORCHECK_unsupported_op( const char*     _file_,
                                const int       _line_,
                                const char*     _func_ );

/*! FUNCTION:  ERRORCHECK_stacktrace()
 *  SYNOPSIS:  Print stacktrace.
 */
void ERRORCHECK_stacktrace();

/*! FUNCTION:  ERRORCHECK_exit()
 *  SYNOPSIS:  Exit program and print stacktrace.
 */
void ERRORCHECK_exit( int exit_flag );

#endif /* _ERROR_HANDLER_H */