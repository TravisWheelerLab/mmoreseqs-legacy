/*******************************************************************************
 *  FILE:      float.c
 *  PURPOSE:   XXX Object
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _XXX_H
#define _XXX_H

/*
 *  FUNCTION:  XXX_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *    RETURN:  Pointer to <buf>
 */
char* XXX_To_String( const XXX   d,
                     char*       buf );

/*
 *  FUNCTION:  XXX_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
int XXX_Compare(  const XXX   a, 
                  const XXX   b );

#endif /* _XXX_H */
