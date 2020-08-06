/*******************************************************************************
 *  FILE:      trace.c
 *  PURPOSE:   TRACE Object
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
#include "structs.h"
#include "utilities.h"
#include "objects.h"

/* header */
#include "trace.h"

/*
 *  FUNCTION:  TRACE_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
inline
char* TRACE_To_String(  const TRACE   	d,
                     	char*       	buf )
{
   sprintf( buf, "{ [%d] (%d,%d) }", d.st, d.i, d.j );
   return buf;
}

/*
 *  FUNCTION:  TRACE_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
inline
int TRACE_Compare( const TRACE     a, 
                   const TRACE     b )
{
   if ( (a.i - b.i) != 0 ) {
      return (a.i - b.i);
   }

   if ( (a.j - b.j) != 0 ) {
      return (a.j - b.j);
   }

   return (a.st - b.st);
}