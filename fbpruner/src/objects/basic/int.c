/*******************************************************************************
 *  FILE:      float.c
 *  PURPOSE:   INT Object
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/* local imports */
#include "../structs.h"
#include "../../utilities/utilities.h"
#include "../objects.h"

/* header */
#include "int.h"

/*
 *  FUNCTION:  INT_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
inline
char* INT_To_String( const INT   d,
                     char*       buf )
{
   sprintf( buf, "%d", d );
   return buf;
}

/*
 *  FUNCTION:  INT_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
inline
int INT_Compare(  const INT   a, 
                  const INT   b )
{
   return (a - b);
}