/*******************************************************************************
 *  FILE:      int.c
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
#include "../../utilities/_utilities.h"
#include "../_objects.h"

/* header */
#include "int.h"

/*! FUNCTION:  INT_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
inline
INT
INT_Create( const INT   data )
{
   return data;
}

/*! FUNCTION:  INT_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
inline
INT
INT_Destroy( INT   data )
{
   return data;
}

/*! FUNCTION:  INT_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
inline
char* 
INT_To_String( const INT   d,
               char*       buf )
{
   sprintf( buf, "%d", d );
   return buf;
}

/*! FUNCTION:  INT_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
inline
int 
INT_Compare(   const INT   a, 
               const INT   b )
{
   return (a - b);
}

/*! FUNCTION:  INT_Swap()
 *  SYNOPSIS:  Swap values of <a> and <b>.
 */
inline
void
INT_Swap(   INT*    a,
            INT*    b )
{
   INT tmp;
   tmp   = *a;
   *a    = *b;
   *b    = tmp;
}