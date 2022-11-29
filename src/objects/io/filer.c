/*******************************************************************************
 *  - FILE:      filer.c
 *  - DESC:    FILER Class. Helps with reading from files.
 *  NOTES:
 *    - Should operations be probably be safer?  Optimization shouldn't be too important with reading/writing.
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
#include "filer.h"


/*!  FUNCTION:  FILER_Create()
 *   SYNOPSIS:  Create <filer>, allocate memory and return pointer.
 */
FILER*
FILER_Create(STR filename,
             STR mode) {
  FILER* filer;
  filer = ERROR_malloc(sizeof(FILER));

  filer->fp = NULL;
  filer->filename = NULL;
  filer->filename = STR_Create(filename);
  /* TODO: need safety check if mode is valid */
  filer->mode = NULL;
  filer->mode = STR_Create(mode);
  filer->is_eof = false;
  filer->is_open = false;

  return filer;
}

/*!  FUNCTION:  FILER_Destroy()
 *   SYNOPSIS:  Destroy <filer>, free memory and return NULL pointer.
 */
FILER*
FILER_Destroy(FILER* filer) {
  /* if already destroyed, do nothing */
  if (filer == NULL) {
    return filer;
  }

  /* close file now if still open */
  if (filer->is_open == true) {
    FILER_Close(filer);
  }
  filer->filename = STR_Destroy(filer->filename);
  filer->mode = STR_Destroy(filer->mode);
  filer = ERROR_free(filer);

  return filer;
}

/*!  FUNCTION:  FILER_Open()
 *   SYNOPSIS:  Open <filer> file pointer.
 */
STATUS_FLAG
FILER_Open(FILER* filer) {
  /* if filename or mode is NULL, opening with throw error */
  if (filer->filename == NULL || filer->mode == NULL) {
    fprintf(stderr, "ERROR: Failed to open file: ( Filename: %s, Mode: %s )\n", filer->filename, filer->mode);
    ERRORCHECK_exit(STATUS_FAILURE);
  }
  /* check if file is already open or if it is stdout (in which case it is always left open) */
  if (filer->is_open || filer->fp == stdout) {
    return STATUS_FAILURE; /* need flag for double-open/double-close */
  }

  filer->fp = fopen(filer->filename, filer->mode);
  if (filer->fp == NULL) {
    fprintf(stderr, "NULL FILEPOINTER IN INDEX READER_Open\n");
    exit(1);
  }

  filer->is_open = true;
  filer->is_eof = false;

  return STATUS_SUCCESS;
}

/*  FUNCTION:  FILER_Close()
 *  SYNOPSIS:  Close <filer> file pointer.
 */
STATUS_FLAG
FILER_Close(FILER* filer) {
  /* special case: don't close standard output file pointers */
  if (FILER_Is_StandardOutput(filer) == false) {
    ERROR_fclose(filer->fp);
    filer->fp = NULL;
    filer->is_open = false;
  }

  return STATUS_SUCCESS;
}

/*!  FUNCTION:  FILER_Is_EndOfFile()
 *   SYNOPSIS:  Check if <filer> is at END_OF_FILE.
 */

/*!  FUNCTION:  FILER_Is_Open()
 *   SYNOPSIS:  Check if <filer> is at open or closed.
 */

/*!  FUNCTION:  FILER_Is_StandardOutput()
 *   SYNOPSIS:  Check if <filer> is to standard output (either stdout or stderr).
 */
bool FILER_Is_StandardOutput(FILER* filer) {
  bool is_stdout, is_stderr;
  is_stdout = (filer->fp == stdout);
  is_stderr = (filer->fp == stderr);
  return is_stdout || is_stderr;
}
