/*******************************************************************************
 *  FILE:      char.c
 *  PURPOSE:   CHAR Object
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _CHAR_H
#define _CHAR_H

/*
 *  FUNCTION:  CHAR_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
char* CHAR_To_String( 	const CHAR   	d,
                     	char*       	buf );

/*
 *  FUNCTION:  CHAR_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
int CHAR_Compare(  	const CHAR   a, 
                  	const CHAR   b );

#endif /* _CHAR_H */
