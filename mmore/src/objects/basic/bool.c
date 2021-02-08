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

/*! FUNCTION:  BOOL_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
inline
char* 
BOOL_To_String( const BOOL   data,
               char*       buf )
{
   sprintf( buf, "%s", data ? "T" : "F" );
   return buf;
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
BOOL_Compare(   const BOOL   a, 
               const BOOL   b )
{
   return (a - b);
}