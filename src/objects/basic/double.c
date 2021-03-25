/*******************************************************************************
 *  FILE:      double.c
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
#include "../../utilities/_utilities.h"
#include "../_objects.h"

/* header */
#include "_basic.h"
#include "double.h"

/*! FUNCTION:  DBL_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
inline
DBL
DBL_Create( const DBL   data )
{
   return data;
}

/*! FUNCTION:  DBL_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
inline
DBL
DBL_Destroy( DBL   data )
{
   return data;
}

/*! FUNCTION:  DBL_Clear()
 *  SYNOPSIS:  Clear <data>.  If pointer data, sets to null. Otherwise, do nothing. 
 */
inline
DBL
DBL_Clear( DBL   data )
{
   return data;
}

/*! FUNCTION:  DBL_ToString()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
inline
char* 
DBL_ToString(  const DBL   data,
               char*       buf )
{
   sprintf( buf, "%.3f", data );
   return buf;
}

/*! FUNCTION:  DBL_FromString()
 *  SYNOPSIS:  Extracts data from string.
 *
 *    RETURN:  Pointer to <buf>
 */
inline
DBL
DBL_FromString(  char*   str )
{
   DBL data;
   data = atof( str );
   return data; 
}

/*! FUNCTION:  DBL_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
inline
int 
DBL_Compare(   const DBL   a, 
               const DBL   b )
{
   return (a - b);
}

/*! FUNCTION:  DBL_CompareTo()
 *  SYNOPSIS:  Generic compare. Casts then compares <a> and <b>.
 *    RETURN:  POS if (a > b), 
 *             0 if equal, 
 *             NEG if (a < b)
 */
inline
int 
DBL_CompareTo(    const void*   a, 
                  const void*   b )
{
   DBL* x = (DBL*)a;
   DBL* y = (DBL*)b;

   return DBL_Compare( *x, *y );
}
