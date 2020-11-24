/*******************************************************************************
 *  FILE:      string.c
 *  PURPOSE:   STRING utility functions.
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
#include "structs.h"
#include "../utilities/utilities.h"
#include "objects.h"

/* header */
#include "mystring.h"


char*
STRING_Create( char*          dest,
               const char*    str ) 
{
   
}

/* Get the number of characters in a string (including \0) */
int 
STRING_Len( const char*  str ) 
{
   int i;
   for (i = 0; str[i] != '\0'; i++) {}
   return i;
}

/* Return concatenation of two input strings as new string */
char* 
STRING_Concat( const char*   str1,
               const char*   str2 )
{
   char*    str;
   size_t   len;

   len = strlen(str1) + strlen(str2) + 1;
   str = (char*) malloc( len * sizeof(char) );

   if ( str == NULL ) {
      fprintf(stderr, "ERROR: malloc failed in STRING_Concat.\n");
      exit(EXIT_FAILURE);
   }

   strcpy(str, str1);
   strcat(str, str2);

   return str;
}

/* replace old character with new character in string (NOTE: Cannot be a string literal) */
void 
STRING_Replace(   char*          str,
                  const char     ch_old,
                  const char     ch_new )
{
   char* p = NULL;
   int i = 0;
   for ( p = str; *p != '\0'; p++ ) {
      if ( *p == ch_old ) *p = ch_new;
   }
}

/* terminates string at given char */
void STRING_Terminate_At(  char*       str,
                           const char  ch )
{
   char* p = NULL;
   int i = 0;
   for ( p = str; *p != '\0'; p++ ) {
      if ( *p == ch ) *p = '\0';
   }
}

/* compares first n chars of string 1 and string 2 */ 
int 
STRING_BeginsWith(   const char*     str1,
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
int 
STRING_EndsWith(  const char*    str1,
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

/* returns filename of full filepath */
char* 
STRING_GetFileFromPath( const char* in_filepath ) 
{
   char*    filepath    = strdup( in_filepath );
   char*    token       = NULL;
   char*    prv_token   = NULL;
   char*    delim       = "/";

   /* get last non-null token in string while breaking on backslash */

   while( ( token = strtok( filepath, delim ) ), token != NULL  ) {
      prv_token = token;
   }

   ERROR_free(filepath);
   return prv_token;
}

/* returns location of search_string in string. */
int 
STRING_Find(   const char* str,
               const char* search )
{
   return 0;
}

/* returns true or false if strings are equal */
int 
STRING_Equal(  const char* str1,
               const char* str2 )
{
   if ( str1 == NULL ) return false;
   if ( str2 == NULL ) return false;

   if ( strcmp( str1, str2 ) == 0 ) {
      return true;
   } else {
      return false;
   }
}