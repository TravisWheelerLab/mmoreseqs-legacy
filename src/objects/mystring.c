/*******************************************************************************
 *  FILE:      string.c
 *  PURPOSE:   STRING utility functions
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
#include "../../utility.h"

/* header */
#include "mystring.h"

/* Get the number of characters in a string (including \0) */
int STRING_Len( const char*  str ) 
{
   int i;
   for (i = 0; str[i] != '\0'; i++) {}
   return i;
}

/* compares first n chars of string 1 and string 2 */ 
int STRING_BeginsWith( const char*     str1,
                       const char*     str2,
                       const size_t    n )
{
   int len1 = strlen( str1 );
   int len2 = strlen( str2 );

   if ( len1 < n || len2 < n ) return -1;

   for ( int i = 0; i < n; i++ ) {
      if ( str1[i] != str2[i] ) {
         return (str1[i] - str2[i]);
      }
   }

   return 0;
}

/* compares last n chars of string 1 and string 2 */
int STRING_EndsWith( const char*    str1,
                     const char*    str2,
                     const size_t   n )
{
   int len1 = strlen( str1 );
   int len2 = strlen( str2 );

   if ( len1 < n || len2 < n ) return -1;

   // printf( "COMPARE>: %s vs %s\n", str1, str2 );
   // printf( "COMPARE: %s vs %s\n", &(str1[len1-(n-1)]), &(str2[len2-(n-1)]) );

   return strncmp( &(str1[len1-(n-1)]), &(str2[len2-(n-1)]), n-1 );
}