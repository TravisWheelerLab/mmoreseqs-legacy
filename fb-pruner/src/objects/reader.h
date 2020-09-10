/*******************************************************************************
 *  FILE:      reader.c
 *  PURPOSE:   READER Class.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _READER_H
#define _READER_H

/*
 *  FUNCTION:  READER_Create()
 *  SYNOPSIS:
 */
READER* READER_Create( STR filename );

/*
 *  FUNCTION:  READER_Destroy()
 *  SYNOPSIS:
 */
void* READER_Destroy( READER* reader );

/*  FUNCTION:  READER_Open()
 *  SYNOPSIS:
 */
int READER_Open( READER* reader );

/*  FUNCTION:  READER_Close()
 *  SYNOPSIS:
 */
int READER_Close( READER* reader );

/*  FUNCTION:  READER_Rewind()
 *  SYNOPSIS:  Set file pointer to beginning of file.
 */
int READER_Rewind( READER* reader );

/*  FUNCTION:  READER_JumpTo()
 *  SYNOPSIS:  Jump to nth byte of into file pointer.
 */
int READER_JumpTo(   READER*     reader, 
                     long int    offset );

/*  FUNCTION:  READER_GetLine()
 *  SYNOPSIS:
 */
char* READER_GetLine( READER* reader );

/*  FUNCTION:  READER_Split()
 *  SYNOPSIS:  Split buffer into 
 */
VECTOR_CHAR* READER_Split( READER* reader );

#endif /* _READER_H */
