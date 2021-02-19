/*******************************************************************************
 *  FILE:      buffer_read.h
 *  PURPOSE:   BUFFER Class. 
 *             Functions for readint from buffer.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:     
 *******************************************************************************/

#ifndef _BUFFER_READ_H
#define _BUFFER_READ_H

/*!  FUNCTION:  BUFFER_GetLine()
 *   SYNOPSIS:  Get current line from <buffer>.
 *   RETURN:    Pointer to head of <line> buffer, NULL if is at END_OF_FILE.
 */
STR
BUFFER_GetLine(   BUFFER*     buffer );

/*!  FUNCTION:  BUFFER_NextLine()
 *   SYNOPSIS:  Get next line from file pointer <fp> and store in <buffer>.
 *   RETURN:    Pointer to head of <line> buffer, NULL if is at END_OF_FILE.
 */
STR
BUFFER_NextLine(  BUFFER*     buffer,
                  FILE*       fp );

/*!  FUNCTION:  BUFFER_RemoveNewline()
 *   SYNOPSIS:  Remove newline from end of <line> buffer.
 */
STATUS_FLAG
BUFFER_RemoveNewline(  BUFFER*     buffer );

/*!  FUNCTION:  BUFFER_SetDelimiter()
 *   SYNOPSIS:  Split <buffer>'s <delimiter> to be used for splitting lines.
 */
STATUS_FLAG
BUFFER_SetDelimiter(    BUFFER*     buffer,
                        STR         new_delimiter );

/*!  FUNCTION:  BUFFER_SplitLine()
 *   SYNOPSIS:  Split <buffer>'s <line> by <delim> and load into <field_offsets>.
 *   RETURN:    Returns number of line splits / fields.
 */
size_t
BUFFER_SplitLine(    BUFFER*     buffer );

/*!  FUNCTION:  BUFFER_SplitLineOn()
 *   SYNOPSIS:  Split <buffer>'s <line> by <delim> and load into <field_offsets>.
 *              Uses temporary <delimiter> supplied by user.
 *   RETURN:    Returns number of line splits / fields.
 */
size_t
BUFFER_SplitLineOn(  BUFFER*     buffer,
                     STR         delimiter );

/*!  FUNCTION:  BUFFER_GetField()
 *   SYNOPSIS:  Get <i>th field in <line> buffer. 
 *              Caller must have called BUFFER_SplitLine().   
 *   RETURN:    Pointer to head of null-terminated string <field>.
 */
STR
BUFFER_GetField(  BUFFER*     buffer,
                  size_t      i );

/*!  FUNCTION:  BUFFER_NextField()
 *   SYNOPSIS:  Get Next Field in <line> buffer.
 *              Caller must have called BUFFER_SplitLine().   
 *   RETURN:    Pointer to head of buffer. If at end of fields, returns NULL.
 */
STR
BUFFER_NextField(    BUFFER*     buffer );

#endif /* _BUFFER_READ_H */