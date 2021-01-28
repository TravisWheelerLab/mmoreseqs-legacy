/*******************************************************************************
 *  FILE:      bound.h
 *  PURPOSE:   BOUND Object
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _BOUND_H
#define _BOUND_H

/*! FUNCTION:  BOUND_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
char* BOUND_To_String(  const BOUND    d,
                        char*          buf );

/*! FUNCTION:  BOUND_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
int BOUND_Compare(   const BOUND   a, 
                  	const BOUND   b );

#endif /* _BOUND_H */
