/*******************************************************************************
 *  FILE:      reader.c
 *  PURPOSE:   READER Class.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:     
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

#endif /* _DEBUGGER_H */