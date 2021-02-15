/*******************************************************************************
 *  FILE:      char.c
 *  PURPOSE:   CHAR Object
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
#include "../../utilities/_utilities.h"
#include "../_objects.h"

/* header */
#include "_basic.h"
#include "char.h"

/*! FUNCTION:  CHAR_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
inline
CHAR
CHAR_Create( const CHAR   data )
{
   return data;
}

/*! FUNCTION:  CHAR_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
inline
CHAR
CHAR_Destroy( CHAR   data )
{
   return data;
}

/*! FUNCTION:  CHAR_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
inline
char* 
CHAR_To_String( 	const CHAR   	d,
                  char*       	buf )
{
   sprintf( buf, "%d", d );
   return buf;
}

/*! FUNCTION:  CHAR_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
inline
int 
CHAR_Compare(  const CHAR   a, 
               const CHAR   b )
{
   return (a - b);
}

/*! FUNCTION:  CHAR_CompareTo()
 *  SYNOPSIS:  Generic compare. Casts then compares <a> and <b>.
 *    RETURN:  POS if (a > b), 
 *             0 if equal, 
 *             NEG if (a < b)
 */
inline
int 
CHAR_CompareTo(    const void*   a, 
                  const void*   b )
{
   CHAR* x = (CHAR*)a;
   CHAR* y = (CHAR*)b;

   return CHAR_Compare( *x, *y );
}
