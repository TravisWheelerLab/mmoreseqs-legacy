/*******************************************************************************
 *  FILE:      bool.c
 *  PURPOSE:   BOOL Object
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
#include "bool.h"

/*! FUNCTION:  BOOL_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
inline
BOOL
BOOL_Create( const BOOL   data )
{
   return data;
}

/*! FUNCTION:  BOOL_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
inline
BOOL
BOOL_Destroy( BOOL   data )
{
   return data;
}

/*! FUNCTION:  BOOL_Clear()
 *  SYNOPSIS:  Clear <data>.  If pointer data, sets to null. Otherwise, do nothing. 
 */
inline
BOOL
BOOL_Clear( BOOL   data )
{
   return data;
}

/*! FUNCTION:  BOOL_ToString()
 *  SYNOPSIS:  Create a string representation of <data>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
inline
char* 
BOOL_ToString( const BOOL     data,
               char*          buf )
{
   sprintf( buf, "%s", data ? "T" : "F" );
   return buf;
}

/*! FUNCTION:  BOOL_FromString()
 *  SYNOPSIS:  Extracts data from string.
 *
 *    RETURN:  Pointer to <buf>
 */
inline
BOOL
BOOL_FromString(  char*   str )
{
   BOOL data;
   data = atoi( str );
   return data; 
}

/*! FUNCTION:  BOOL_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
inline
int 
BOOL_Compare(     const BOOL   a, 
                  const BOOL   b )
{
   return (a - b);
}

/*! FUNCTION:  BOOL_CompareTo()
 *  SYNOPSIS:  Generic compare. Casts then compares <a> and <b>.
 *    RETURN:  POS if (a > b), 
 *             0 if equal, 
 *             NEG if (a < b)
 */
inline
int 
BOOL_CompareTo(    const void*   a, 
                  const void*   b )
{
   BOOL* x = (BOOL*)a;
   BOOL* y = (BOOL*)b;

   return BOOL_Compare( *x, *y );
}
