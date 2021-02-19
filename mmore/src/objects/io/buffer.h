/*******************************************************************************
 *  FILE:      buffer.h
 *  PURPOSE:   BUFFER Class. Helps with reading from files.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:     
 *******************************************************************************/

#ifndef _BUFFER_H
#define _BUFFER_H

/*!  FUNCTION:  BUFFER_Create()
 *   SYNOPSIS:  Create <buffer>, allocate memory and return pointer.  
 */
BUFFER*
BUFFER_Create();

/*!  FUNCTION:  BUFFER_Destroy()
 *   SYNOPSIS:  Destroy <buffer>, free memory and return NULL pointer.  
 */
BUFFER*
BUFFER_Destroy( BUFFER*  buffer );

/*!  FUNCTION:  BUFFER_SetDelimiter()
 *   SYNOPSIS:  Split <buffer>'s <delimiter> to be used for splitting lines.
 */
STATUS_FLAG
BUFFER_SetDelimiter(    BUFFER*     buffer,
                        STR         new_delimiter );

/*!  FUNCTION:  BUFFER_Move()
 *   SYNOPSIS:  Move <i> signed positions from current position in <buffer>.
 *              WARNING: Does not do bounds checking.
 */
STR
BUFFER_Move(   BUFFER*     buffer,
               int         i );

/*!  FUNCTION:  BUFFER_Move()
 *   SYNOPSIS:  Checks whether buffer is empty.
 */
bool
BUFFER_IsEmpty(   BUFFER*     buffer );

#endif /* _BUFFER_H */