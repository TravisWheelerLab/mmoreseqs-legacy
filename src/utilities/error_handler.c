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

void ERRORCHECK_Handler( const int    error_code,
                    const char*  file,
                    const int    line,
                    const char*  func,
                    const char*  note )
{
    switch ( error_code )
    {
    case ERROR_UNKNOWN:
        fprintf( stderr, "%s\n", "ERROR: Unknown error." );
        break;
    case ERROR_MALLOC:
        fprintf( stderr, "%s\n", "ERROR: Malloc error." );
        break;
    case ERROR_REALLOC:
        fprintf(stderr, "%s\n", "ERROR: Realloc error." );
        break;
    default:
        fprintf( stderr, "%s\n", "ERROR: Error occurred with invalid error code." );
    }
    if (file != NULL)
        fprintf(stderr, "ERROR occurred in FILE: \"%s\", at LINE: %d, in FUNC: \"%s\".\n)", file, line, func);
    fprintf(stderr, "[ CODE = %d ]\n", error_code);
    if (note != NULL)
        fprintf(stderr, "NOTE: %s\n", note);

    fprintf(stderr, "Program terminated.\n");
    exit(EXIT_FAILURE);
}

/* check for null pointer */
void ERRORCHECK_NullPtr( const void*   data,
                         int*          error_code,
                         const char*   note )
{

}

/* open file and check for null file pointer */ 
FILE* ERRORCHECK_fopen( char*       filename,
                        const char* permission )
{
    FILE* fp = fopen( filename, permission );
    if (fp == NULL) {
        fprintf(stderr, "ERROR: Unable to open file '%s' for '%s'.\n", filename, permission);
    }
}

/* malloc and check for null pointer */
void* ERRORCHECK_malloc()
{

}

/* report error location */
void ERRORCHECK_location( const char* _file_,
                          const int   _line_,
                          const char* _func_  )
{
    fprintf(stderr, "ERROR occurred in FILE: %s, LINE: %d, FUNC: %s.\n", _file_, _line_, _func_);
}

/* error if called function is currently unsupported */
void ERRORCHECK_unsupported_op( )
{

}