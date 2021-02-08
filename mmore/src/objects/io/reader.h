/*******************************************************************************
 *  FILE:      reader.c
 *  PURPOSE:   READER Class.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _READER_H
#define _READER_H

/*!  FUNCTION:  READER_Create()
 *   SYNOPSIS:  Create <reader>, allocate memory and return pointer.  
 */
READER* 
READER_Create( STR filename );

/*!  FUNCTION:  READER_Destroy()
 *   SYNOPSIS:  Destroy <reader>, free memory and return NULL pointer.
 */
READER* 
READER_Destroy( READER* reader );

/*!  FUNCTION:  READER_Open()
 *   SYNOPSIS:  Open <reader> file pointer.
 */
STATUS_FLAG 
READER_Open( READER* reader );

/*!  FUNCTION:  READER_Close()
 *   SYNOPSIS:  Close <reader> file pointer.
 */
STATUS_FLAG 
READER_Close( READER* reader );

/*!  FUNCTION:  READER_Rewind()
 *   SYNOPSIS:  Sets <reader> file pointer to beginning of file.
 */
STATUS_FLAG 
READER_Rewind( READER* reader );

/*!  FUNCTION:  READER_JumpTo()
 *   SYNOPSIS:  Jump <reader> file pointer to <offset>th byte of file.
 */
STATUS_FLAG 
READER_JumpTo(    READER*     reader, 
                  long int    offset );

/*!  FUNCTION:  READER_GetLine()
 *   SYNOPSIS:  Get current line from <reader>.
 *   RETURN:    Pointer to head of <line> buffer, NULL if is at END_OF_FILE.
 */
STR
READER_GetLine(   READER*  reader );

/*!  FUNCTION:  READER_NextLine()
 *   SYNOPSIS:  Get next line from file <fp>. 
 */
STR
READER_NextLine(  READER*     reader );

/*!  FUNCTION:  READER_SetDelimiter()
 *   SYNOPSIS:  Split <buffer>'s <delimiter> to be used for splitting lines.
 */
STATUS_FLAG
READER_SetDelimiter(    READER*     reader,
                        STR         new_delimiter );

/*!  FUNCTION:  READER_SplitLine()
 *   SYNOPSIS:  Split <reader>'s <line> by <delim> and load into <field_offsets>.
 *   RETURN:    Returns number of line splits / fields.
 */
size_t
READER_SplitLine(    READER*     reader );

/*!  FUNCTION:  READER_GetField()
 *   SYNOPSIS:  Get <i>th field in <line> buffer. 
 *              Caller must have called READER_SplitLine().   
 *   RETURN:    Pointer to head of null-terminated <field> in <line>.
 */
STR
READER_GetField(  READER*     reader,
                  size_t      i );

/*!  FUNCTION:  READER_NextField()
 *   SYNOPSIS:  Get Next Field in <line> buffer.
 *              Caller must have called READER_SplitLine().   
 *   RETURN:    Pointer to head of buffer. If at end of fields, returns NULL.
 */
STR
READER_NextField(    READER*     reader );

/*!  FUNCTION:  READER_Is_EndOfFile()
 *   SYNOPSIS:  Check if <reader> is at END_OF_FILE.
 */
bool
READER_Is_EndOfFile(    READER*     reader );

#endif /* _READER_H */
