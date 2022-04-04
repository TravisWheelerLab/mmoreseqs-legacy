/*******************************************************************************
 *  - FILE:      buffer.c
 *  - DESC:    BUFFER Class. Helps with reading from files.
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

/*!  FUNCTION:  BUFFER_Create()
 *   SYNOPSIS:  Create <buffer>, allocate memory and return pointer.
 */
BUFFER*
BUFFER_Create() {
  BUFFER* buffer;
  buffer = ERROR_malloc(sizeof(BUFFER));

  buffer->line = VECTOR_CHAR_Create();
  buffer->split_line = VECTOR_CHAR_Create();
  buffer->line_len = -1;

  buffer->field = NULL;
  buffer->field_offsets = VECTOR_INT_Create();
  buffer->delim = STR_Create("\t");
  buffer->field_idx = -1;

  return buffer;
}

/*!  FUNCTION:  BUFFER_Destroy()
 *   SYNOPSIS:  Destroy <buffer>, free memory and return NULL pointer.
 */
BUFFER*
BUFFER_Destroy(BUFFER* buffer) {
  int line_max = 2048;
  if (buffer == NULL) {
    return buffer;
  }

  buffer->line_max = line_max;
  buffer->line = VECTOR_CHAR_Destroy(buffer->line);
  buffer->split_line = VECTOR_CHAR_Destroy(buffer->split_line);
  buffer->field_offsets = VECTOR_INT_Destroy(buffer->field_offsets);
  buffer->delim = STR_Destroy(buffer->delim);
  buffer = ERROR_free(buffer);

  return buffer;
}

/*!  FUNCTION:  BUFFER_Empty()
 *   SYNOPSIS:  Empty <buffer>
 */
void BUFFER_Empty(BUFFER* buffer) {
  VECTOR_CHAR_Reuse(buffer->line);
  VECTOR_CHAR_Set(buffer->line, 0, NULL_CHAR);
}

/*!  FUNCTION:  BUFFER_SetDelimiter()
 *   SYNOPSIS:  Split <buffer>'s <delimiter> to be used for splitting lines.
 */
STATUS_FLAG
BUFFER_SetDelimiter(BUFFER* buffer,
                    STR new_delimiter) {
  buffer->delim = STR_Set(buffer->delim, new_delimiter);
  return STATUS_SUCCESS;
}

/*!  FUNCTION:  BUFFER_Move()
 *   SYNOPSIS:  Move <i> signed positions from current position in <buffer>.
 *              WARNING: Does not do bounds checking.
 */
STR BUFFER_Move(BUFFER* buffer,
                int i) {
  buffer->buffer_ptr += i;
  return buffer->buffer_ptr;
}

/*!  FUNCTION:  BUFFER_IsEmpty()
 *   SYNOPSIS:  Checks whether buffer is empty.
 */
bool BUFFER_IsEmpty(BUFFER* buffer) {
  bool is_empty;
  is_empty = (buffer->line_len <= 0);
  return is_empty;
}

/*!  FUNCTION:  BUFFER_GetLength()
 *   SYNOPSIS:  Get occupied length of buffer.
 */
size_t
BUFFER_GetLength(BUFFER* buffer) {
  return VECTOR_CHAR_GetSize(buffer->line);
}

/*!  FUNCTION:  BUFFER_GetMaxLength()
 *   SYNOPSIS:  Get max length of buffer.
 */
size_t
BUFFER_GetMaxLength(BUFFER* buffer) {
  return VECTOR_CHAR_GetSize(buffer->line);
}
