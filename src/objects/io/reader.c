/*******************************************************************************
 *  - FILE:      reader.c
 *  - DESC:    READER Class. Helps with reading from files.
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
#include "reader.h"

/*!  FUNCTION:  READER_Create()
 *   SYNOPSIS:  Create <reader>, allocate memory and return pointer.
 */
READER*
READER_Create(STR filename) {
  READER* reader;
  reader = ERROR_malloc(sizeof(READER));

  reader->fp = NULL;
  reader->filename = STR_Create(filename);
  reader->buffer = BUFFER_Create();
  reader->file_size = 0;
  reader->is_eof = false;
  reader->is_open = false;
  reader->line_count = 0;
  reader->lines_read = 0;

  return reader;
}

/*!  FUNCTION:  READER_Destroy()
 *   SYNOPSIS:  Destroy <reader>, free memory and return NULL pointer.
 */
READER*
READER_Destroy(READER* reader) {
  if (reader == NULL) {
    return reader;
  }

  if (reader->fp != NULL) {
    READER_Close(reader);
  }
  reader->filename = STR_Destroy(reader->filename);
  reader->buffer = BUFFER_Destroy(reader->buffer);
  reader = ERROR_free(reader);

  return reader;
}

/*!  FUNCTION:  READER_Open()
 *   SYNOPSIS:  Open <reader> file pointer.
 */
STATUS_FLAG
READER_Open(READER* reader) {
  reader->fp = fopen(reader->filename, "r");

  if (reader->fp == NULL) {
    fprintf(stderr, "NULL FILEPOINTER IN INDEX READER_Open\n");
    exit(1);
  }

  reader->is_open = true;
  reader->is_eof = false;

  return STATUS_SUCCESS;
}

/*  FUNCTION:  READER_Close()
 *  SYNOPSIS:
 */
STATUS_FLAG
READER_Close(READER* reader) {
  ERROR_fclose(reader->fp);
  reader->fp = NULL;
  reader->is_open = false;

  return STATUS_SUCCESS;
}

/*!  FUNCTION:  READER_Rewind()
 *   SYNOPSIS:  Sets <reader> file pointer to beginning of file.
 */
STATUS_FLAG
READER_Rewind(READER* reader) {
  STATUS_FLAG status;
  status = fseek(reader->fp, 0, SEEK_SET);
  return status;
}

/*!  FUNCTION:  READER_JumpTo()
 *   SYNOPSIS:  Jump <reader> file pointer to <offset>th byte of file.
 */
STATUS_FLAG
READER_JumpTo(READER* reader,
              long int offset) {
  STATUS_FLAG status;
#if SAFE || TRUE
  {
    if (reader->is_open == false) {
      fprintf(stderr, "ERROR: Trying to JumpTo() before opening file.\n");
      ERRORCHECK_exit(EXIT_FAILURE);
    }
  }
#endif
  status = fseek(reader->fp, offset, SEEK_SET);
  return status;
}

/*!  FUNCTION:  READER_GetLine()
 *   SYNOPSIS:  Get current line from <reader>.
 *   RETURN:    Pointer to head of <line> buffer, NULL if is at END_OF_FILE.
 */
STR READER_GetLine(READER* reader) {
  STR line;
  line = BUFFER_GetLine(reader->buffer);
  return line;
}

/*!  FUNCTION:  READER_NextLine()
 *   SYNOPSIS:  Get next line from file <fp>.
 */
STR READER_NextLine(READER* reader) {
  STR line;
  /* read line */
  BUFFER_NextLine(reader->buffer, reader->fp);
  BUFFER_RemoveNewline(reader->buffer);
  line = BUFFER_GetLine(reader->buffer);
  /* update stats */
  reader->pos_prv = reader->pos_cur;
  reader->pos_cur = ftell(reader->fp);
  reader->lines_read += 1;
  reader->line_count += 1;
  reader->is_eof = (line == NULL);

  return line;
}

/*!  FUNCTION:  READER_SetDelimiter()
 *   SYNOPSIS:  Split <buffer>'s <delimiter> to be used for splitting lines.
 */
STATUS_FLAG
READER_SetDelimiter(READER* reader,
                    STR new_delimiter) {
  STATUS_FLAG status;
  status = BUFFER_SetDelimiter(reader->buffer, new_delimiter);
  return status;
}

/*!  FUNCTION:  READER_SplitLine()
 *   SYNOPSIS:  Split <reader>'s <line> by <delim> and load into <field_offsets>.
 *   RETURN:    Returns number of line splits / fields.
 */
size_t
READER_SplitLine(READER* reader) {
  size_t num_fields;
  num_fields = BUFFER_SplitLine(reader->buffer);
  return num_fields;
}

/*!  FUNCTION:  READER_GetField()
 *   SYNOPSIS:  Get <i>th field in <line> buffer.
 *              Caller must have called READER_SplitLine().
 *   RETURN:    Pointer to head of null-terminated <field> in <line>.
 */
STR READER_GetField(READER* reader,
                    size_t i) {
  STR field;
  field = BUFFER_GetField(reader->buffer, i);
  return field;
}

/*!  FUNCTION:  READER_NextField()
 *   SYNOPSIS:  Get Next Field in <line> buffer.
 *              Caller must have called READER_SplitLine().
 *   RETURN:    Pointer to head of buffer. If at end of fields, returns NULL.
 */
STR READER_NextField(READER* reader) {
  STR field;
  field = BUFFER_NextField(reader->buffer);
  return field;
}

/*!  FUNCTION:  READER_Is_EndOfFile()
 *   SYNOPSIS:  Check if <reader> is at END_OF_FILE.
 */
bool READER_Is_EndOfFile(READER* reader) {
  return reader->is_eof;
}
