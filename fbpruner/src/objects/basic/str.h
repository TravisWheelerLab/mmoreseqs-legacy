/*******************************************************************************
 *  FILE:      str.c
 *  PURPOSE:   STR Object ( wraps *char )
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _STR_H
#define _STR_H

/*
 *  FUNCTION:  FLT_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
char* 
STR_To_String( const STR   d,
               char*       buf );

/*
 *  FUNCTION:  FLT_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
int 
STR_Compare(  	const STR   a, 
               const STR   b );

#endif /* _STR_H */
