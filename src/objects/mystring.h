/*******************************************************************************
 *  FILE:      string.h
 *  PURPOSE:   STRING utility functions
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _STRING_H
#define _STRING_H

/* Get the number of characters in a string (including \0) */
int STRING_Len( const char*  str );

/* Return concatenation of two input strings as new string */
char* STRING_Concat( const char*   str1,
                     const char*   str2 );

/* compares first n chars of string 1 and string 2 */ 
int STRING_BeginsWith( const char*     str1,
                       const char*     str2,
                       const size_t    n );

/* compares last n chars of string 1 and string 2 */
int STRING_EndsWith( const char*    str1,
                     const char*    str2,
                     const size_t   n );

#endif /* _STRING_H */