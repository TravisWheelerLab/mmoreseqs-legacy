/*******************************************************************************
 *  FILE:      buffer_write.c
 *  PURPOSE:   BUFFER Class. 
 *             Functions for writing to buffer.
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
#include <sys/types.h>
#include <sys/stat.h>

/* local imports */
#include "../structs.h"
#include "../../utilities/_utilities.h"
#include "../_objects.h"

/* header */
#include "_io.h"
#include "buffer_write.h"

/*!  FUNCTION:  BUFFER_Write()
 *   SYNOPSIS:  Write STRING <str> to BUFFER <buffer> (only if it does not overflow buffer).
 *              Return true if buffer is full.
 */
bool
BUFFER_Write(  BUFFER*     buffer,
               STR         str )
{
   int N, N_buf, N_str, N_max; 
   
   N_buf = BUFFER_GetLength( buffer );
   N_str = STR_GetLength( str );
   N_max = BUFFER_GetMaxLength( buffer );
   /* verify buffer is large enough to store new string */
   if ( N_buf + N_str > N_max ) {
      return true;
   }
   /* add string to buffer */
   VECTOR_CHAR_Append( buffer->line, str, N_str );
   return false;
}

/*!  FUNCTION:  BUFFER_WriteLine()
 *   SYNOPSIS:  Write STRING <str> to BUFFER <buffer>, with a newline character at the end.
 *              Return true if buffer is full.
 */
bool
BUFFER_WriteLine(    BUFFER*     buffer,
                     STR         str )
{
   int N, N_buf, N_str, N_max; 
   
   N_buf = BUFFER_GetLength( buffer );
   N_str = STR_GetLength( str );
   N_max = BUFFER_GetMaxLength( buffer );
   /* verify buffer is large enough to store new string */
   if ( N_buf + N_str + 1 > N_max ) {
      return true;
   }
   /* add string to buffer */
   VECTOR_CHAR_Append( buffer->line, str, N_str );
   VECTOR_CHAR_Push( buffer->line, '\n' );
   return false;
}

