/*******************************************************************************
 *  FILE:      filer.c
 *  PURPOSE:   FILER Class. 
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _FILER_H
#define _FILER_H

/*!  FUNCTION:  FILER_Create()
 *   SYNOPSIS:  Create <filer>, allocate memory and return pointer.  
 */
FILER* 
FILER_Create(  STR filename,
               STR mode );

/*!  FUNCTION:  FILER_Destroy()
 *   SYNOPSIS:  Destroy <filer>, free memory and return NULL pointer.
 */
FILER* 
FILER_Destroy( FILER* filer );

/*!  FUNCTION:  FILER_SetFile()
 *   SYNOPSIS:  Open <filer> file pointer.
 */
STATUS_FLAG
FILER_SetFile( FILER*   filer,
               STR      filename,
               STR      mode );

/*!  FUNCTION:  FILER_Open()
 *   SYNOPSIS:  Open <filer> file pointer.
 */
STATUS_FLAG 
FILER_Open( FILER* filer );

/*!  FUNCTION:  FILER_Close()
 *   SYNOPSIS:  Close <filer> file pointer.
 */
STATUS_FLAG 
FILER_Close( FILER* filer );

/*!  FUNCTION:  FILER_Rewind()
 *   SYNOPSIS:  Sets <filer> file pointer to beginning of file.
 */
STATUS_FLAG 
FILER_Rewind( FILER* filer );

/*!  FUNCTION:  FILER_JumpTo()
 *   SYNOPSIS:  Jump <filer> file pointer to <offset>th byte of file.
 */
STATUS_FLAG 
FILER_JumpTo(    FILER*     filer, 
                  long int    offset );

/*!  FUNCTION:  FILER_NextLine()
 *   SYNOPSIS:  Get next line from file <fp>. 
 */
STR
FILER_NextLine(  FILER*     filer );

/*!  FUNCTION:  FILER_Is_EndOfFile()
 *   SYNOPSIS:  Check if <filer> is at END_OF_FILE.
 */
bool
FILER_Is_EndOfFile(    FILER*     filer );

/*!  FUNCTION:  FILER_Is_Open()
 *   SYNOPSIS:  Check if <filer> is at open or closed.
 */
bool
FILER_Is_Open(    FILER*     filer );

/*!  FUNCTION:  FILER_Is_StandardOutput()
 *   SYNOPSIS:  Check if <filer> is to standard output.
 */
bool
FILER_Is_StandardOutput(    FILER*     filer );

#endif /* _FILER_H */
