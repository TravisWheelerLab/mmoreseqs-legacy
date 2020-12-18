/*******************************************************************************
 *  FILE:      float.c
 *  PURPOSE:   DBL Object
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _DBL_H
#define _DBL_H

/*
 *  FUNCTION:  DBL_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
char* DBL_To_String( const DBL   d,
                     char*       buf );

/*
 *  FUNCTION:  DBL_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
int DBL_Compare(  const DBL   a, 
                  const DBL   b );

#endif /* _DBL_H */
