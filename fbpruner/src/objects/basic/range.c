/*******************************************************************************
 *  FILE:      trace.c
 *  PURPOSE:   RANGE Object
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
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
#include "range.h"

/*
 *  FUNCTION:  RANGE_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
inline
char* RANGE_To_String(  const RANGE   	d,
                     	char*       	buf )
{
   sprintf( buf, "(%d,%d)", d.beg, d.end );
   return buf;
}

/*
 *  FUNCTION:  RANGE_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
inline
int RANGE_Compare( const RANGE     a, 
                   const RANGE     b )
{
   if ( (a.beg - b.beg) != 0 ) {
      return (a.beg - b.beg);
   }

   return (a.end - b.end);
}