/*******************************************************************************
 *  FILE:      xstr.c
 *  PURPOSE:   XSTR object.
 *             String object which tracks its size and supports more safe string functions.
 *             Not optimal for string building.
 *
 *  AUTHOR:    Dave Rich
 *  NOTES:
 *    - WIP.       
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
#include "xstr.h"

/*! FUNCTION:  X_STRING_Create()
 *  SYNOPSIS:  Allocates data to hold <xstr> with length of <src>.
 */
XSTR* 
XSTR_Create( char*   src )
{
   XSTR* str;
   str         = ERROR_malloc( sizeof(XSTR) );
   str->data   = strdup( src );
   str->N      = strlen( src );

   return str;
}

/*! FUNCTION:  XSTR_Destroy()
 *  SYNOPSIS:  Destroys <str> data and frees memory.
 */
XSTR* 
XSTR_Destroy( XSTR* str )
{
   ERROR_free( str->data );
   ERROR_free( str );
   return NULL;
}

/*! FUNCTION:  XSTR_Copy()
 *  SYNOPSIS:  Creates a deep copy of <src>, stores it in <dest> and returns it.
 */
XSTR* 
XSTR_Copy(  XSTR*    src,
            XSTR*    dest )
{
   if ( dest == NULL ) {
      dest = XSTR_Create( src->data );
      return dest;
   }

   dest->data  = ERROR_free( dest->data );
   dest->data  = strdup( src->data );
   dest->N     = src->N;

   return dest;
}

/*! FUNCTION:  XSTR_SetSize()
 *  SYNOPSIS:  Allocates <str> to store a string of length <L>.
 *             
 */
XSTR* 
XSTR_SetSize( XSTR*      str,
               size_t     L )
{
   str->data   = ERROR_realloc( str->data, sizeof(char) * (L+1) );
   str->N      = L;

   return str;
}

/*! FUNCTION:  X_STRING_Set()
 *  SYNOPSIS:  Set <str> to hold <src>.
 */
XSTR* 
XSTR_Set( XSTR*   str,
          char*   src )
{
   if ( src == NULL ) return NULL;

   
   str->data   = ERROR_free( str->data );
   str->data   = strdup( src );
   str->N      = strlen( src );

   return str;
}

/*! FUNCTION:  XSTR_Concat()
 *  SYNOPSIS:  Copies <str_1> and <str_2> (in order) into new <str>.
 *             Returns <str>.
 */
XSTR* 
XSTR_Concat( XSTR*   str_1,
             XSTR*   str_2 )
{
   XSTR* str      = ERROR_malloc( sizeof(XSTR) );

   str->data      = ERROR_malloc( sizeof(char) * str_1->N );
   str->data[0]   = '\0'; 
   strcat( str->data, str_1->data );
   strcat( str->data, str_2->data );

   return str;
}

/*! FUNCTION:  XSTR_Append()
 *  SYNOPSIS:  Appends <str_add> <str> and <str_2> (in order) into new <str>.
 *             Returns <str>.
 */
XSTR* 
XSTR_Append( XSTR*   str,
             XSTR*   str_add )
{
   XSTR_SetSize( str, str->N + str_add->N );
   str->data[0]   = '\0'; 
   strcat( str->data, str_add->data );

   return str;
}