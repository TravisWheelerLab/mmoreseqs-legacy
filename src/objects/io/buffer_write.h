/*******************************************************************************
 *  - FILE:      buffer_write.h
 *  - DESC:    BUFFER Class.
 *             Functions for writing to buffer.
 *******************************************************************************/

#ifndef _BUFFER_WRITE_H
#define _BUFFER_WRITE_H

/*!  FUNCTION:  BUFFER_Write()
 *   SYNOPSIS:  Write STRING <str> to BUFFER <buffer>.
 *              Return true if buffer is full.
 */
bool BUFFER_Write(BUFFER* buffer, STR str);

/*!  FUNCTION:  BUFFER_WriteLine()
 *   SYNOPSIS:  Write STRING <str> to BUFFER <buffer>, with a newline character
 * at the end. Return true if buffer is full.
 */
bool BUFFER_WriteLine(BUFFER* buffer, STR str);

#endif /* _BUFFER_WRITE_H */
