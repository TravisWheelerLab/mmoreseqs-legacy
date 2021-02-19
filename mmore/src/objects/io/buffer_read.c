/*******************************************************************************
 *  FILE:      buffer.c
 *  PURPOSE:   BUFFER Class. Helps with reading from files.
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
#include "buffer.h"

/*!  FUNCTION:  BUFFER_GetLine()
 *   SYNOPSIS:  Get current line from <buffer>.
 *   RETURN:    Pointer to head of <line> buffer, NULL if is at END_OF_FILE.
 */
STR
BUFFER_GetLine(   BUFFER*     buffer )
{
   if ( BUFFER_IsEmpty(buffer) ) return NULL;

   STR line;
   line = buffer->line->data;
   return line;
}

/*!  FUNCTION:  BUFFER_NextLine()
 *   SYNOPSIS:  Get next line from file pointer <fp> and store in <buffer>.
 *   RETURN:    Pointer to head of <line> buffer, NULL if is at END_OF_FILE.
 */
STR
BUFFER_NextLine(  BUFFER*     buffer,
                  FILE*       fp )
{
   /* get next line for <fp> and store string and length in <line> buffer. */
   buffer->line_len     = getline( &(buffer->line->data), &(buffer->line->Nalloc), fp );
   buffer->line->N      = buffer->line_len;
   /* set buffer pointer to head of buffer */
   buffer->buffer_ptr   = buffer->line->data;
   /* if <line> buffer is empty (zero length), then buffer is at the END_OF_FILE */
   bool is_eof = ( buffer->line_len <= 0 );
   if ( is_eof ) {
      STR_SetEmpty( buffer->line->data );
      return NULL;
   }
   return buffer->line->data;
}

/*!  FUNCTION:  BUFFER_RemoveNewline()
 *   SYNOPSIS:  Remove newline from end of <line> buffer.
 */
STATUS_FLAG
BUFFER_RemoveNewline(  BUFFER*     buffer )
{
   bool is_empty;
   bool ends_with_newline;

   is_empty = BUFFER_IsEmpty(buffer);
   /* if buffer is not empty... */
   if (  is_empty == false ) {
      /* and buffer ends in newline character */
      ends_with_newline = VEC_X( buffer->line, buffer->line->N - 1 ) == '\n';
      if ( ends_with_newline == true ) {
         /* then remove it from the buffer */
         VEC_X( buffer->line, buffer->line->N - 1 ) = '\0';
         buffer->line->N -= 1;
      }
   }
}

/*!  FUNCTION:  BUFFER_SplitLine()
 *   SYNOPSIS:  Split <buffer>'s <line> by <delim> and load into <field_offsets>.
 *   RETURN:    Returns number of line splits / fields.
 */
size_t
BUFFER_SplitLine(    BUFFER*     buffer )
{
   size_t num_fields;

   /* reinit field data */
   VECTOR_INT_Reuse( buffer->field_offsets );
   buffer->field_idx = 0;

   /* TODO: copy or just reference? */
   VECTOR_CHAR_Copy( buffer->line, buffer->split_line );
   STR   split_line  = buffer->split_line->data;
   STR   field       = split_line;
   STR   delimiter   = buffer->delim;

   /* loop through <line> and split into fields by delimiter */
   field = strtok( field, delimiter );
   while ( field != NULL ) 
   {
      int offset = field - split_line;
      VECTOR_INT_Pushback( buffer->field_offsets, offset );
      field = strtok( field, delimiter );
   }
   num_fields = VECTOR_INT_GetSize( buffer->field_offsets );
   return num_fields;
}

/*!  FUNCTION:  BUFFER_SplitLineOn()
 *   SYNOPSIS:  Split <buffer>'s <line> by <delim> and load into <field_offsets>.
 *              Uses temporary <delimiter> supplied by user.
 *   RETURN:    Returns number of line splits / fields.
 */
size_t
BUFFER_SplitLineOn(  BUFFER*     buffer,
                     STR         delimiter )
{
   size_t   num_fields;
   STR      tmp_delimiter = delimiter;

   /* swap contents of given delimiter and default delimiter, split string, then swap back */
   STR_Swap( &tmp_delimiter, &buffer->delim );
   num_fields = BUFFER_SplitLine( buffer );
   STR_Swap( &tmp_delimiter, &buffer->delim );

   return num_fields;
}

/*!  FUNCTION:  BUFFER_GetField()
 *   SYNOPSIS:  Get <i>th field in <line> buffer. 
 *              Caller must have called BUFFER_SplitLine().   
 *   RETURN:    Pointer to head of null-terminated <field> in <line>.
 */
STR
BUFFER_GetField(  BUFFER*     buffer,
                  size_t      i )
{
   int offset;
   STR field;
   offset  = VEC_X( buffer->field_offsets, i );
   field   = VEC_XX( buffer->split_line, offset );
   return field;
}

/*!  FUNCTION:  BUFFER_NextField()
 *   SYNOPSIS:  Get Next Field in <line> buffer.
 *              Caller must have called BUFFER_SplitLine().   
 *   RETURN:    Pointer to head of buffer. If at end of fields, returns NULL.
 */
STR
BUFFER_NextField(    BUFFER*     buffer )
{
   STR field;
   /* check if we've reached the end of the fields in line */
   if ( buffer->field_idx >= buffer->field_offsets->N ) {
      return NULL;
   }
   /* if not, get field */
   field = BUFFER_GetField( buffer, buffer->field_idx );
   buffer->field_idx += 1;
   return field;
}

/*!  FUNCTION:  BUFFER_GetFieldLength()
 *   SYNOPSIS:  Get length of <i>th <field> in <buffer>.
 */
size_t
BUFFER_GetFieldLength(  BUFFER*     buffer,
                        size_t      i )
{
   STR      field;
   size_t   len;
   field    = BUFFER_GetField( buffer, i );
   len      = STR_GetLength( field );
   return len;
}