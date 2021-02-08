/*******************************************************************************
 *  FILE:      str.c
 *  PURPOSE:   STR Object.
 *             Mostly a simple wrapper for char*.
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

/* header */
#include "../_objects.h"
#include "_basic.h"
#include "str.h"

/*! FUNCTION:  STR_Create()
 *  SYNOPSIS:  Create <str>, fills with input <str_chrs> and return pointer to <str>.
 */
inline
STR
STR_Create(  const char*   chrs )
{
   /* if string points to null, return null */
   if (chrs == NULL) return NULL;

   STR str;
   str = strdup( chrs );
   return str;
}

/*! FUNCTION:  STR_Destroy()
 *  SYNOPSIS:  Destroys <str>, frees memory and returns NULL pointer.
 */
inline
STR
STR_Destroy( STR   str )
{
   if ( str == NULL ) return str;

   str = ERROR_free( str );

   return str;
}

/*! FUNCTION:  STR_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *    RETURN:  Pointer to <buf>
 */
inline
char* 
STR_To_String( const STR   data,
               char*       buf )
{
   sprintf( buf, "%s", data );

   return buf;
}

/*! FUNCTION:  STR_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
inline
int 
STR_Compare(  const STR   a, 
              const STR   b )
{
   return strcmp( a, b );
}