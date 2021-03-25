/*******************************************************************************
 *  FILE:      writer.c
 *  PURPOSE:   WRITER Class.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _WRITER_H
#define _WRITER_H

/*!  FUNCTION:  WRITER_Create()
 *   SYNOPSIS:  Create <writer>, allocate memory and return pointer.  
 */
WRITER* 
WRITER_Create( STR filename );

/*!  FUNCTION:  WRITER_Destroy()
 *   SYNOPSIS:  Destroy <writer>, free memory and return NULL pointer.
 */
WRITER* 
WRITER_Destroy( WRITER* writer );

/*!  FUNCTION:  WRITER_Open()
 *   SYNOPSIS:  Open <writer> file pointer.
 */
STATUS_FLAG 
WRITER_Open( WRITER* writer );

/*!  FUNCTION:  WRITER_Close()
 *   SYNOPSIS:  Close <writer> file pointer.
 */
STATUS_FLAG 
WRITER_Close( WRITER* writer );

/*!  FUNCTION:  WRITER_Rewind()
 *   SYNOPSIS:  Sets <writer> file pointer to beginning of file.
 */
STATUS_FLAG 
WRITER_Rewind( WRITER* writer );

/*!  FUNCTION:  WRITER_JumpTo()
 *   SYNOPSIS:  Jump <writer> file pointer to <offset>th byte of file.
 */
STATUS_FLAG 
WRITER_JumpTo(    WRITER*     writer, 
                  long int    offset );

#endif /* _WRITER_H */
