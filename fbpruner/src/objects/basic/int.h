/*******************************************************************************
 *  FILE:      float.c
 *  PURPOSE:   INT Object
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _INT_H
#define _INT_H

/*
 *  FUNCTION:  INT_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *    RETURN:  Pointer to <buf>
 */
char* INT_To_String( const INT   d,
                     char*       buf );

/*
 *  FUNCTION:  INT_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
int INT_Compare(  const INT   a, 
                  const INT   b );

#endif /* _INT_H */
