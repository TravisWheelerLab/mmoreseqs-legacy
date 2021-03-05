/*******************************************************************************
 *  FILE:      str.c
 *  PURPOSE:   STR Object.  
 *             Extended functionality.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *    - None known.
 *  NOTES:
 *    - This will contain functions that are beyond those in the other basic types.
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
#include "str_ext.h"

/* null terminator character */
const char null_char = '\0';

/*! FUNCTION:  STR_Create_with_Size()
 *  SYNOPSIS:  Create <str> that can store str of size <L> and return pointer.
 *             No changes made to allocated data.
 */
inline
STR
STR_Create_with_Size(  const size_t    L )
{
   STR str;
   str = ERROR_malloc( sizeof(char) * (L+1) );

   return str;
}

/*! FUNCTION:  STR_Create_Empty_with_Size()
 *  SYNOPSIS:  Create empty <str> with size enough to store string of length <L> (adds room for null character).
 *             Use with caution, as it does not track size of allocated memory.
 */
inline
STR
STR_Create_Empty_with_Size( const size_t  L )
{
   STR str;
   str = STR_Create_with_Size( L );
   STR_SetEmpty( str );

   return str;
}

/*! FUNCTION:  STR_Set()
 *  SYNOPSIS:  Sets <str> to be copy of contents in <str_chrs>.
 *             Destroys <str_old>, creates a string to fills with input <str_new> to replace it, 
 *             and return pointer.  <str_old> should be reassigned to returned pointer.
 */
inline
STR
STR_Set( STR            str_old,
         const char*    chr_new )
{
   if (chr_new == NULL) {
      str_old  = STR_Destroy( str_old );
      return NULL;
   }
   
   STR str_new;
   str_new  = STR_Create( chr_new );
   str_old  = STR_Destroy( str_old );
   return str_new;
}

/*! FUNCTION:  STR_Copy()
 *  SYNOPSIS:  Copies <src> to <dest>.
 *             Destroys <str_old>, reates a string to fills with input <str_new> to replace it, 
 *             and return pointer.  <str_old> should be reassigned to returned pointer.
 */
inline
STR
STR_Copy(   STR            dest,
            const STR      src )
{
   STR str;
   dest  = STR_Destroy( dest );
   str   = STR_Create( src );
   return str;
}

/*! FUNCTION:  STR_Equals()
 *  SYNOPSIS:  Compares data_1 and data_2. Returns true if they are equal, false otherwise.
 */
inline
bool
STR_Equals(    const STR   data_1,
               const STR   data_2 )
{
   int cmp;
   bool is_equal;
   cmp      = STR_Compare( data_1, data_2 );
   is_equal = ( cmp == 0 ); 
   return is_equal;
}

/*! FUNCTION:  STR_Overwrite_Char()
 *  SYNOPSIS:  Overwrites the contents of <ch> to <src>, and returns <dest>.
 */
inline
STR
STR_Overwrite_Char(  STR      dest,
                     char     ch )
{
   dest[0] = ch;
   return dest;
}

/*! FUNCTION:  STR_Overwrite_Str()
 *  SYNOPSIS:  Overwrites the contents of <str> to <dest>, and returns <dest>.
 *  WARNING:   Use with caution, does no memory checks.
 */
inline
STR
STR_Overwrite_Str(   STR      dest,
                     STR      str )
{
   dest = strcpy( dest, str );
   return dest;
}

/*! FUNCTION:  STR_SetEmpty()
 *  SYNOPSIS:  Make <str> an empty string. Sets first character to null.
 *             Use with caution, as it does not track size of allocated memory.
 */
inline
STR
STR_SetEmpty( STR  str )
{
   str[0] = null_char;

   return str;
}

/*! FUNCTION:  STR_Concat()
 *  SYNOPSIS:  Creates new string <str> that contains a concatenation of <str_1> and <str_2>.
 *             Returns <str>.
 */
inline
STR
STR_Concat(    const STR   str_1,
               const STR   str_2 )
{
   if ( str_1 == NULL ) {
      return STR_Create( str_2 );
   }
   if ( str_2 == NULL ) {
      return STR_Create( str_1 );
   }

   STR      str;
   size_t   L;

   L     = STR_GetLength( str_1 ) + STR_GetLength( str_2 );
   str   = STR_Create_Empty_with_Size( L );
   str   = strcat( str, str_1 );
   str   = strcat( str, str_2 );

   return str;
}

/*! FUNCTION:  STR_Append()
 *  SYNOPSIS:  Appends <str_add> to the end of <str>.  
 *             Returns <str>.
 *   WARNING:  Reallocates <str>, so <str> should be reassigned to returned <str> value, 
 *             as input <str> will be freed.
 */
inline
STR
STR_Append(    STR         str,
               const STR   str_add )
{
   STR      str_old;
   STR      str_new;
   
   str_old  = str;
   str_new  = STR_Concat( str_old, str_add );
   str_old  = STR_Destroy( str_old );
   str      = str_new;

   return str;
}

/*! FUNCTION:  STR_GetLength()
 *  SYNOPSIS:  Returns length of <str>.
 */
inline
size_t
STR_GetLength( const STR   str )
{
   if ( str == NULL ) {
      return 0;
   }

   size_t L;
   L = strlen( str );
   return L;
}

/*! FUNCTION:  STR_GetX()
 *  SYNOPSIS:  Gets location of <i>th character in <str>.
 */
inline
STR
STR_GetX(   STR         str,
            const int   i )
{
   STR substr;
   substr = &(str[i]);
   return substr;
}

/*! FUNCTION:  STR_GetChar()
 *  SYNOPSIS:  Gets <i>th character in <str>.
 */
inline
char
STR_GetChar(   const STR   str,
               const int   i )
{
   char ch;
   ch = str[i];
   return ch;
}

/*! FUNCTION:  STR_Find_Char()
 *  SYNOPSIS:  Linear search for first instance of <ch> in <str> and returns <str_loc>, 
 *             which points to location in <str>.
 *             Returns <str_loc>, or NULL if not in <str>.
 */
inline
STR
STR_Find_Char(    STR      str,
                  char     ch )
{
   STR str_loc;
   str_loc = strchr( str, ch );
   return str_loc;
}

/*! FUNCTION:  STR_Find_Substring()
 *  SYNOPSIS:  Linear search for <substr> in <str> and returns <loc>.
 *             Returns <loc>, or -1 if not in <str>.
 */
inline
STR
STR_Find_Substring(  STR   str,
                     STR   substr )
{
   STR str_loc;
   str_loc = strstr( str, substr );
   return str_loc;
}

/*! FUNCTION:  STR_Replace_Char()
 *  SYNOPSIS:  Replace every instance of <old_char> in string and replaces it with <new_char>.
 */
STR
STR_Replace_Char(    STR      str,
                     char     ch_old,
                     char     ch_new )
{
   STR   str_find    = str;

   /* get first instance */
   str_find = STR_Find_Char( str, ch_old );
   /* repeat until no more instances found */
   while ( str_find != NULL ) 
   {
      /* replace <ch_old> with <ch_new> */
      STR_Overwrite_Char( str_find, ch_new );
      /* find next instance of <ch_old> in <str> */
      str_find += 1;
      str_find  = STR_Find_Char( str_find, ch_old );
   } 
   return str;
}

/*! FUNCTION:  STR_Replace_Substring()
 *  SYNOPSIS:  Replace every instance of <old_char> in string and replaces it with <new_char>.
 */
inline
STR
STR_Replace_Substring(     STR     str,
                           STR     substr_old,
                           STR     substr_new )
{
   STR   str_find    = str;

   /* get first instance */
   str_find = STR_Find_Substring( str, substr_old );
   /* repeat until no more instances found */
   while ( str_find != NULL ) 
   {
      /* replace <ch_old> with <ch_new> */
      STR_Overwrite_Str( str_find, substr_new );
      /* find next instance of <ch_old> in <str> */
      str_find += 1;
      str_find  = STR_Find_Substring( str_find, substr_old );
   } 
   return str;
}

/*! FUNCTION:  STR_Compare_Prefix()
 *  SYNOPSIS:  Compare first <L> characters of <a> and <b>.
 *    RETURN:  POS if (a > b), 
 *             0 if equal, 
 *             NEG if (a < b)
 */
int 
STR_Compare_Prefix(  const STR      a, 
                     const STR      b,
                     const size_t   L )
{
   return strncmp( a, b, L );
}

/*! FUNCTION:  STR_Swap()
 *  SYNOPSIS:  Swap strings pointed to by <str_a> and <str_b>.
 */
STATUS_FLAG 
STR_Swap(   STR*      str_a, 
            STR*      str_b )
{
   STR   str_tmp;
   str_tmp  = *str_a;
   *str_a   = *str_b;
   *str_b   = str_tmp;

   return STATUS_SUCCESS;
}

/*! FUNCTION:  STR_ToUpper()
 *  SYNOPSIS:  Convert <str> to all upper case.
 */
STR 
STR_ToUpper(   STR      str )
{
   int i = 0;
   while ( str[i] != NULL_CHAR ) {
      str[i] = toupper(str[i]);
      i++;
   }
   return str;
}

/*! FUNCTION:  STR_ToLower()
 *  SYNOPSIS:  Convert <str> to all lower case.
 */
STR 
STR_ToLower(   STR      str )
{
   int i = 0;
   while ( str[i] != NULL_CHAR ) {
      str[i] = tolower(str[i]);
      i++;
   }
   return str;
}

/*! FUNCTION:  STR_StartsWith()
 *  SYNOPSIS:  Checks if start of <str> and <prefix> are equal up to the length of <prefix>.
 */
bool  
STR_StartsWith(   STR      str,
                  STR      prefix )
{
   int   is_startswith;
   int   S  = STR_GetLength( str );
   int   P  = STR_GetLength( prefix );

   /* if <prefix> is longer than <str>, then <str> cannot contain it */
   if ( P > S ) return -1;

   /* otherwise, compare against the length of <prefix>. */
   is_startswith = strncmp( str, prefix, P );
   return is_startswith;
}

/*! FUNCTION:  STR_EndsWith()
 *  SYNOPSIS:  Checks if end of <str> and <suffix> are equal up to the length of <suffix>.
 */
bool  
STR_EndsWith(  STR      str,
               STR      suffix )
{
   int   is_endswith;
   int   S  = STR_GetLength( str );
   int   P  = STR_GetLength( suffix );

   /* if <prefix> is longer than <str>, then <str> cannot contain it */
   if ( P > S ) return -1;

   /* otherwise, compare against the length of <prefix>. */
   str = &(str[S - P]);
   is_endswith = strncmp( str, suffix, P );
   return is_endswith;
}

