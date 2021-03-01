/*******************************************************************************
 *  FILE:       ptr.c
 *  PURPOSE:    PTR Object.
 *              Void pointer wrapper.
 *
 *  AUTHOR:     Dave Rich
 *  BUG:
 *    - None known.   
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
#include "ptr.h"

/*! FUNCTION:  PTR_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
inline
PTR
PTR_Create( const PTR   data )
{
   return data;
}

/*! FUNCTION:  PTR_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
inline
PTR
PTR_Destroy( PTR   data )
{
   return data;
}

/*! FUNCTION:  PTR_Empty()
 *  SYNOPSIS:  Clear <data>.  If pointer data, sets to null. Otherwise, do nothing. 
 */
inline
PTR
PTR_Clear( PTR   data )
{
   return data;
}

/*! FUNCTION:  PTR_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *    RETURN:  POS if (a > b), 
 *             0 if equal, 
 *             NEG if (a < b)
 */
inline
int 
PTR_Compare(   const PTR   a, 
               const PTR   b )
{
   return (a - b);
}

/*! FUNCTION:  PTR_CompareTo()
 *  SYNOPSIS:  Generic compare. Casts then compares <a> and <b>.
 *    RETURN:  POS if (a > b), 
 *             0 if equal, 
 *             NEG if (a < b)
 */
inline
int 
PTR_CompareTo(    const void*   a, 
                  const void*   b )
{
   PTR* x = (PTR*)a;
   PTR* y = (PTR*)b;

   return PTR_Compare( *x, *y );
}

/*! FUNCTION:  PTR_Swap()
 *  SYNOPSIS:  Swap values of <a> and <b>.
 */
inline
void
PTR_Swap(   PTR*    a,
            PTR*    b )
{
   PTR tmp;
   tmp   = *a;
   *a    = *b;
   *b    = tmp;
}

/*! FUNCTION:  PTR_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a preallocated char* buffer <buf> of length <buf_size>.
 *    RETURN:  Pointer to <buf>.
 */
inline
char* 
PTR_To_String(    const PTR   data,
                  char*       buf )
{
   sprintf( buf, "%p", data );
   return buf;
}
