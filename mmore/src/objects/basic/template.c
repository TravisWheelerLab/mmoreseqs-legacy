/*******************************************************************************
 *  FILE:       template.c
 *  PURPOSE:    XXX Object
 *
 *  AUTHOR:     Dave Rich
 *  BUG:
 *    - None known.   
 *  NOTE: 		
 *    - Generated using TEMPLATE via "scripts/builder-helper/build_basic_wrapper_files.sh" 
 *    - These wrappers contain some superfluous methods.  This just fulfills an API for simpler vectors.   
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
#include "template.h"

/*! FUNCTION:  XXX_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
inline
XXX
XXX_Create( const XXX   data )
{
   return data;
}

/*! FUNCTION:  XXX_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
inline
XXX
XXX_Destroy( XXX   data )
{
   return data;
}

/*! FUNCTION:  XXX_Clear()
 *  SYNOPSIS:  Clear <data>.  If pointer data, sets to null. Otherwise, do nothing. 
 */
inline
XXX
XXX_Clear( XXX   data )
{
   return data;
}

/*! FUNCTION:  XXX_Swap()
 *  SYNOPSIS:  Swap values of <a> and <b>.
 */
inline
void
XXX_Swap(   XXX*    a,
            XXX*    b )
{
   XXX tmp;
   tmp   = *a;
   *a    = *b;
   *b    = tmp;
}

/*! FUNCTION:  XXX_ToString()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a preallocated char* buffer <buf> of length <buf_size>.
 *    RETURN:  Pointer to <buf>.
 */
inline
char* 
XXX_ToString(    const XXX   data,
                  char*       buf )
{
   sprintf( buf, "%d", data );
   return buf;
}

/*! FUNCTION:  XXX_FromString()
 *  SYNOPSIS:  Extracts data from string.
 *
 *    RETURN:  Pointer to <buf>
 */
inline
XXX
XXX_FromString(  char*   str )
{
   XXX data;
   data = atoi( str );
   return data; 
}

/*! FUNCTION:  XXX_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *    RETURN:  POSITIVE if (a > b), 
 *             ZERO if equal, 
 *             NEGATIVE if (a < b)
 */
inline
int 
XXX_Compare(   const XXX   a, 
               const XXX   b )
{
   return (a - b);
}

/*! FUNCTION:  XXX_CompareTo()
 *  SYNOPSIS:  Generic compare. Casts then compares <a> and <b>.
 *    RETURN:  POSITIVE if (a > b), 
 *             ZERO if equal, 
 *             NEGATIVE if (a < b)
 */
inline
int 
XXX_CompareTo(    const void*   a, 
                  const void*   b )
{
   XXX* x = (XXX*)a;
   XXX* y = (XXX*)b;

   return XXX_Compare( *x, *y );
}

/*! FUNCTION:  XXX_Equals()
 *  SYNOPSIS:  Checks if <a> and <b> are equal.
 *    RETURN:  TRUE if equal, FALSE otherwise
 */
inline
bool 
XXX_Equals(    const XXX   a, 
               const XXX   b )
{
   bool equals;
   equals = ( XXX_Compare(a,b) == 0 );
   return equals;
}
