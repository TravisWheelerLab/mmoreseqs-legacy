/*******************************************************************************
 *  FILE:      buffer_write.h
 *  PURPOSE:   BUFFER Class. 
 *             Functions for writing to buffer.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:     
 *******************************************************************************/

#ifndef _BUFFER_WRITE_H
#define _BUFFER_WRITE_H

/*!  FUNCTION:  BUFFER_Write()
 *   SYNOPSIS:  Write STRING <str> to BUFFER <buffer>.
 *              Return true if buffer is full.
 */
bool
BUFFER_Write(   BUFFER*     buffer );

/*!  FUNCTION:  BUFFER_WriteLine()
 *   SYNOPSIS:  Write STRING <str> to BUFFER <buffer>, with a newline character at the end.
 *              Return true if buffer is full.
 */
bool
BUFFER_WriteLine(   BUFFER*     buffer );

#endif /* _BUFFER_WRITE_H */