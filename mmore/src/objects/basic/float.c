/*******************************************************************************
 *  FILE:      float.c
 *  PURPOSE:   FLT Object
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
#include "float.h"

/*! FUNCTION:  FLT_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
inline
FLT
FLT_Create( const FLT   data )
{
   return data;
}

/*! FUNCTION:  FLT_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
inline
FLT
FLT_Destroy( FLT   data )
{
   return data;
}

/*! FUNCTION:  FLT_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
inline
char* 
FLT_To_String( const FLT   d,
               char*       buf )
{
   sprintf( buf, "%.3f", d );
   return buf;
}

/*! FUNCTION:  FLT_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
inline
int 
FLT_Compare(   const FLT   a, 
               const FLT   b )
{
   return (a - b);
}

/*! FUNCTION:  FLT_CompareTo()
 *  SYNOPSIS:  Generic compare. Casts then compares <a> and <b>.
 *    RETURN:  POS if (a > b), 
 *             0 if equal, 
 *             NEG if (a < b)
 */
inline
int 
FLT_CompareTo(    const void*   a, 
                  const void*   b )
{
   FLT* x = (FLT*)a;
   FLT* y = (FLT*)b;

   return FLT_Compare( *x, *y );
}
