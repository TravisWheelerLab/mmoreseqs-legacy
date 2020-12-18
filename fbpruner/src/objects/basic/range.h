/******************************************************************************
 *  FILE:      trace.c
 *  PURPOSE:   RANGE Object
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _RANGE_H
#define _RANGE_H

/*
 *  FUNCTION:  RANGE_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
char* RANGE_To_String(  const RANGE   	d,
                     	char*       	buf );

/*
 *  FUNCTION:  RANGE_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
int RANGE_Compare( const RANGE  a, 
                   const RANGE  b );

#endif /* _RANGE_H */