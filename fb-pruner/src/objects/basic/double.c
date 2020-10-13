/*******************************************************************************
 *  FILE:      float.c
 *  PURPOSE:   DBL Object
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
#include "double.h"

/*
 *  FUNCTION:  DBL_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
inline
char* DBL_To_String( const DBL   d,
                     char*       buf )
{
   sprintf( buf, "%.3f", d );
   return buf;
}

/*
 *  FUNCTION:  DBL_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
inline
int DBL_Compare(  const DBL   a, 
                  const DBL   b )
{
   return (a - b);
}
