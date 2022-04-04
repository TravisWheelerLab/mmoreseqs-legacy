/*******************************************************************************
 *  - FILE:      error_handler.c
 *  - DESC:    Functions for handling error codes.
 *  NOTES:
 *    -  All exits should come through the error handler, hopefully with useful
 *       decriptions if something went wrong.
 *    -  For now, debugging is easier with stacktrace created by -fsanitize=address.
 *       So exits with error codes are sent to stderr and
 *       program allowed to continue if in -BUILD=DEBUG.
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
#include "../objects/_objects.h"

/* header */
#include "_utilities.h"
#include "error_handler.h"

/*!  FUNCTION:  ERRORCHECK_handler()
 *   SYNOPSIS:  Handles errors.
 */
void ERRORCHECK_handler(const ERROR_FLAG error_code,
                        const char* _file_,
                        const int _line_,
                        const char* _func_,
                        const char* err_msg) {
  /* if no message, output default message according to error code */
  if (err_msg == NULL) {
    fprintf(stderr, "ERROR (%d): ", error_code);
    switch (error_code) {
      case ERROR_UNKNOWN: {
        fprintf(stderr, "%s\n", "Unknown error.");
      } break;
      case ERROR_MALLOC: {
        fprintf(stderr, "%s\n", "Malloc error.");
      } break;
      case ERROR_REALLOC: {
        fprintf(stderr, "%s\n", "Realloc error.");
      } break;
      case ERROR_FILE_IO: {
        fprintf(stderr, "%s\n", "File I/O error.");
      }
      default: {
        fprintf(stderr, "%s\n", "Error occurred with invalid error code.");
      }
    }
  }
  /* report custom message */
  else {
    fprintf(stderr, "# ERROR: %s\n", err_msg);
  }
  /* report location of error */
  if (_file_ != NULL && _func_ != NULL) {
    fprintf(stderr, "# ERROR occurred ==> in - FILE: \"%s\", at LINE: %d, in FUNC: \"%s\".\n", _file_, _line_, _func_);
  }
  /* terminate program */
  fprintf(stderr, "# Program terminated.\n");
  ERRORCHECK_exit(EXIT_FAILURE);
}

/*!  FUNCTION:  ERRORCHECK_fopen()
 *   SYNOPSIS:  Opens file <fp> and handles potential errors.
 */
inline FILE*
ERRORCHECK_fopen(const char* filename,
                 const char* permission,
                 const char* _file_,
                 const int _line_,
                 const char* _func_) {
  FILE* fp = fopen(filename, permission);
  if (fp == NULL) {
    ERRORCHECK_handler(ERROR_FILE_IO, _file_, _line_, _func_, filename);
  }
}

/*!  FUNCTION:  ERRORCHECK_fclose()
 *   SYNOPSIS:  Closes file <fp> and handles potential errors.
 */
inline FILE*
ERRORCHECK_fclose(FILE* fp,
                  const char* _file_,
                  const int _line_,
                  const char* _func_) {
  /* close file in not null */
  if (fp != NULL) {
    fclose(fp);
  }
  /* else, attempting to close a NULL file pointer */
  else {
    printf("WARNING: Attempted to close a NULL pointer.");
  }

  return NULL;
}

/*!  FUNCTION:  ERRORCHECK_alloc()
 *   SYNOPSIS:  Reallocates memory <ptr> or size <size> and handles potential errors.
 */
inline void*
ERRORCHECK_alloc(void* ptr,
                 const size_t size,
                 const char* _file_,
                 const int _line_,
                 const char* _func_) {
  /* NOTE: realloc behaves like malloc when {ptr} is NULL. */
  ptr = realloc(ptr, size);
  /* if {ptr} is NULL, then we have a memory error */
  if (ptr == NULL) {
    printf("ERROR_SIZE: %ld\n", size);
    ERRORCHECK_handler(ERROR_MALLOC, _file_, _line_, _func_, NULL);
  }
  return ptr;
}

/*!  FUNCTION:  ERRORCHECK_alloc()
 *   SYNOPSIS:  Allocates memory <ptr> or size <size> and handles potential errors.
 *              Returns new pointer to allocated memory.
 */
inline void*
ERRORCHECK_malloc(const size_t size,
                  const char* _file_,
                  const int _line_,
                  const char* _func_) {
  void* ptr;
  ptr = malloc(size);
  /* if {ptr} is NULL, then we have a memory error */
  if (ptr == NULL) {
    printf("ERROR_SIZE: %ld\n", size);
    ERRORCHECK_handler(ERROR_MALLOC, _file_, _line_, _func_, NULL);
  }
  return ptr;
}

/*!  FUNCTION:  ERRORCHECK_alloc()
 *   SYNOPSIS:  Reallocates memory <ptr> or size <size> and handles potential errors.
 *              Returns new pointer to allocated memory.
 */
inline void*
ERRORCHECK_realloc(void* ptr,
                   const size_t size,
                   const char* _file_,
                   const int _line_,
                   const char* _func_) {
  ptr = realloc(ptr, size);
  /* if {ptr} is NULL, then we have a memory error */
  if (ptr == NULL) {
    printf("ERROR_SIZE: %ld\n", size);
    ERRORCHECK_handler(ERROR_REALLOC, _file_, _line_, _func_, NULL);
  }
  return ptr;
}

/*!  FUNCTION:  ERRORCHECK_free()
 *   SYNOPSIS:  Frees memory <ptr> and handles potential errors.
 *              Returns NULL pointer.
 */
inline void*
ERRORCHECK_free(void* ptr,
                const char* _file_,
                const int _line_,
                const char* _func_) {
  if (ptr != NULL) {
    free(ptr);
  }

  return NULL;
}

/*!  FUNCTION:  ERRORCHECK_boundscheck()
 *   SYNOPSIS:  Checks that <idx> access is with array size <max>.
 */
void* ERRORCHECK_boundscheck(int idx,
                             int max,
                             const int size,
                             const char* _file_,
                             const int _line_,
                             const char* _func_) {
  if (idx > max) {
    ERRORCHECK_handler(ERROR_OUT_OF_BOUNDS, _file_, _line_, _func_, NULL);
  }
}

/*!  FUNCTION:  ERRORCHECK_print_location()
 *   SYNOPSIS:  Prints location in source code.
 */
void ERRORCHECK_print_location(FILE* fp,
                               const char* _file_,
                               const int _line_,
                               const char* _func_) {
  fprintf(fp, "{ - FILE: %s, LINE: %d, FUNCTION: %s }\n",
          _file_, _line_, _func_);
}

/*!  FUNCTION:  ERRORCHECK_unsupported_op()
 *   SYNOPSIS:  Handles error if operation is not supported.
 */
void ERRORCHECK_unsupported_op(const char* _file_,
                               const int _line_,
                               const char* _func_) {
}

/*!  FUNCTION:  ERRORCHECK_memcheck()
 *   SYNOPSIS:  Checks if memory at position has proper value.
 */
void ERRORCHECK_memcheck(int row,
                         int col,
                         float mat,
                         float ins,
                         float del) {
  printf("#> ERROR: Memory at position (%d,%d) not cleared. Value = ( %9.4f %9.4f %9.4f )\n",
         row, col, mat, ins, del);
}

/*!  FUNCTION:  ERRORCHECK_stacktrace()
 *   SYNOPSIS:  Performs stack trace in source code.
 */
void ERRORCHECK_stacktrace() {
  void* array[10];
  size_t size;
  char** strings;
  size_t i;

  size = backtrace(array, 10);
  strings = backtrace_symbols(array, size);

  fprintf(stderr, "Obtained %ld stack frames.\n", size);

  for (i = 0; i < size; i++) {
    fprintf(stderr, "%s\n", strings[i]);
  }
  free(strings);
}

/*!  FUNCTION:  ERRORCHECK_exit()
 *   SYNOPSIS:  Exits program is <exit_flag> and reports errors.
 */
void ERRORCHECK_exit(int exit_flag) {
  if (exit_flag == EXIT_SUCCESS) {
    fprintf(stderr, "# PROGRAM terminated successfully.\n");
    exit(exit_flag);
  } else {
    // ERRORCHECK_stacktrace();
    fprintf(stderr, "ERROR: EXIT REACHED with code: %d\n", exit_flag);
    // exit(exit_flag);

    /* force address sanitizer to do stacktrace */
    int test_arr[1];
    test_arr[2] = 42;
  }
}
