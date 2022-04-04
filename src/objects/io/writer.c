/*******************************************************************************
 *  - FILE:      writer.c
 *  - DESC:    WRITER Class. Helps with reading from files.
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
#include "writer.h"

/*!  FUNCTION:  WRITER_Create()
 *   SYNOPSIS:  Create <writer>, allocate memory and return pointer.
 */
WRITER*
WRITER_Create(STR filename) {
  WRITER* writer;
  writer = ERROR_malloc(sizeof(WRITER));

  writer->fp = NULL;
  writer->filename = STR_Create(filename);
  writer->buffer = BUFFER_Create();
  writer->file_size = 0;
  writer->is_eof = false;
  writer->is_open = false;
  writer->line_count = 0;
  writer->lines_written = 0;

  return writer;
}

/*!  FUNCTION:  WRITER_Destroy()
 *   SYNOPSIS:  Destroy <writer>, free memory and return NULL pointer.
 */
WRITER*
WRITER_Destroy(WRITER* writer) {
  if (writer == NULL) {
    return writer;
  }

  if (writer->fp != NULL) {
    WRITER_Close(writer);
  }
  writer->filename = STR_Destroy(writer->filename);
  writer->buffer = BUFFER_Destroy(writer->buffer);
  writer = ERROR_free(writer);

  return writer;
}

/*!  FUNCTION:  WRITER_Open()
 *   SYNOPSIS:  Open <writer> file pointer.
 */
STATUS_FLAG
WRITER_Open(WRITER* writer) {
  writer->fp = ERROR_fopen(writer->filename, "r");
  writer->is_open = true;
  writer->is_eof = false;

  return STATUS_SUCCESS;
}

/*  FUNCTION:  WRITER_Close()
 *  SYNOPSIS:
 */
STATUS_FLAG
WRITER_Close(WRITER* writer) {
  ERROR_fclose(writer->fp);
  writer->fp = NULL;
  writer->is_open = false;

  return STATUS_SUCCESS;
}

/*!  FUNCTION:  WRITER_Rewind()
 *   SYNOPSIS:  Sets <writer> file pointer to beginning of file.
 */
STATUS_FLAG
WRITER_Rewind(WRITER* writer) {
  STATUS_FLAG status;
  status = fseek(writer->fp, 0, SEEK_SET);
  return status;
}

/*!  FUNCTION:  WRITER_JumpTo()
 *   SYNOPSIS:  Jump <writer> file pointer to <offset>th byte of file.
 */
STATUS_FLAG
WRITER_JumpTo(WRITER* writer,
              long int offset) {
  STATUS_FLAG status;
#if SAFE || TRUE
  {
    if (writer->is_open == false) {
      fprintf(stderr, "ERROR: Trying to JumpTo() before opening file.\n");
      ERRORCHECK_exit(EXIT_FAILURE);
    }
  }
#endif
  status = fseek(writer->fp, offset, SEEK_SET);
  return status;
}

/*!  FUNCTION:  WRITER_Flush()
 *   SYNOPSIS:  Jump <writer> file pointer to <offset>th byte of file.
 */
STATUS_FLAG
WRITER_Flush(WRITER* writer) {
  int N_buffer;
  N_buffer = BUFFER_GetLength(writer->buffer);

  BUFFER_Empty(writer->buffer);
}
