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

/*! FUNCTION:  DBL_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
inline
char* 
DBL_To_String( const DBL   data,
               char*       buf )
{
   sprintf( buf, "%.3f", data );
   return buf;
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
