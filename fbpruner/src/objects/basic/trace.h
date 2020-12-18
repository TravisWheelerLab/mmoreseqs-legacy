/******************************************************************************
 *  FILE:      trace.c
 *  PURPOSE:   TRACE Object
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _TRACE_H
#define _TRACE_H

/*
 *  FUNCTION:  TRACE_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
char* TRACE_To_String(  const TRACE   	d,
                     	char*       	buf );

/*
 *  FUNCTION:  TRACE_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
int TRACE_Compare( const TRACE  a, 
                   const TRACE  b );

#endif /* _TRACE_H */