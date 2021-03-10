/*******************************************************************************
 *  FILE:      str_ext.c
 *  PURPOSE:   STR Object.  
 *             Extended functionality.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _STR_EXT_H
#define _STR_EXT_H

/*! FUNCTION:  STR_Create_with_Size()
 *  SYNOPSIS:  Create <str> that can store str of size <L> and return pointer.
 *             No changes made to allocated data.
 */
STR
STR_Create_with_Size(  const size_t    L );

/*! FUNCTION:  STR_Create_Empty_with_Size()
 *  SYNOPSIS:  Create empty <str> with size enough to store string of length <L> (adds room for null character).
 *             Use with caution, as it does not track size of allocated memory.
 */
STR
STR_Create_Empty_with_Size( const size_t  L );

/*! FUNCTION:  STR_Set()
 *  SYNOPSIS:  Sets <str_old> to <chr_new>.
 *             Destroys <str_old>, creates a string to fills with input <str_new> to replace it, 
 *             and return pointer.  <str_old> should be reassigned to returned pointer.
 */
STR
STR_Set( STR            str_old,
         const char*    chr_new );

/*! FUNCTION:  STR_Copy()
 *  SYNOPSIS:  Copies <src> to <dest>.
 *             Destroys <str_old>, reates a string to fills with input <str_new> to replace it, 
 *             and return pointer.  <str_old> should be reassigned to returned pointer.
 */
STR
STR_Copy(   STR            dest,
            const STR      src );

/*! FUNCTION:  STR_Equals()
 *  SYNOPSIS:  Compares data_1 and data_2. Returns true if they are equal, false otherwise.
 */
bool
STR_Equals(    const STR   data_1,
               const STR   data_2 );

/*! FUNCTION:  STR_Overwrite_Char()
 *  SYNOPSIS:  Overwrites the contents of <ch> to <src>, and returns <dest>.
 */
STR
STR_Overwrite_Char(  STR      dest,
                     char     ch );

/*! FUNCTION:  STR_Overwrite_Str()
 *  SYNOPSIS:  Overwrites the contents of <str> to <dest>, and returns <dest>.
 *  WARNING:   Use with caution, does no memory checks.
 */
STR
STR_Overwrite_Str(   STR      dest,
                     STR      str );

/*! FUNCTION:  STR_SetEmpty()
 *  SYNOPSIS:  Make <str> an empty string. Sets first character to null.
 *             Use with caution, as it does not track size of allocated memory.
 */
STR
STR_SetEmpty( STR  str );

/*! FUNCTION:  STR_Concat()
 *  SYNOPSIS:  Creates new string <str> that contains a concatenation of <str_1> and <str_2>.
 *             Returns <str>.
 */

STR
STR_Concat(    const STR   str_1,
               const STR   str_2 );

/*! FUNCTION:  STR_Append()
 *  SYNOPSIS:  Appends <str_add> to the end of <str>.  
 *             Returns <str>.
 *   WARNING:  Reallocates <str>, so <str> should be reassigned to returned <str> value, 
 *             as input <str> may have been freed.
 */
STR
STR_Append(    STR         str,
               const STR   str_add );

/*! FUNCTION:  STR_GetLength()
 *  SYNOPSIS:  Returns length of <str>.
 */
size_t
STR_GetLength(   const STR   str );

/*! FUNCTION:  STR_GetX()
 *  SYNOPSIS:  Gets location of <i>th character in <str>
 */
STR
STR_GetX(   STR         str,
            const int   i );

/*! FUNCTION:  STR_GetChar()
 *  SYNOPSIS:  Gets <i>th character in <str>.
 */
char
STR_GetChar(  const STR   str,
               const int   i );

/*! FUNCTION:  STR_Find_Char()
 *  SYNOPSIS:  Linear search for first instance of <ch> in <str> and returns <str_loc>, 
 *             which points to location in <str>.
 *  RETURN:    Returns <str_loc>, or NULL if not in <str>.
 */
STR
STR_Find_Char(    STR      str,
                  char     ch );

/*! FUNCTION:  STR_Find_Substring()
 *  SYNOPSIS:  Linear search for <substr> in <str> and returns <loc>.
 *  RETURN:    Returns <str_loc>, or NULL if not in <str>.
 */
STR
STR_Find_Substring(  STR   str,
                     STR   substr );

/*! FUNCTION:  STR_Replace_Char()
 *  SYNOPSIS:  Replace every instance of <old_char> in string and replaces it with <new_char>.
 */
STR
STR_Replace_Char(    STR            str,
                     const char     ch_old,
                     const char     ch_new );

/*! FUNCTION:  STR_Replace_Substring()
 *  SYNOPSIS:  Replace every instance of <old_char> in string and replaces it with <new_char>.
 */
STR
STR_Replace_Substring(     STR     str,
                           STR     substr_old,
                           STR     substr_new );

/*! FUNCTION:  STR_ComparePrefix()
 *  SYNOPSIS:  Compare first <L> characters of <a> and <b>.
 *    RETURN:  POS if (a > b), 
 *             0 if equal, 
 *             NEG if (a < b)
 */
int 
STR_ComparePrefix(  const STR      a, 
                     const STR      b,
                     const size_t   L );

/*! FUNCTION:  STR_Swap()
 *  SYNOPSIS:  Swap strings pointed to by <str_a> and <str_b>.
 */
STATUS_FLAG 
STR_Swap(   STR*      str_a, 
            STR*      str_b );

/*! FUNCTION:  STR_ToUpper()
 *  SYNOPSIS:  Convert <str> to all upper case in-place.
 */
STR 
STR_ToUpper(   STR      str );

/*! FUNCTION:  STR_ToLower()
 *  SYNOPSIS:  Convert <str> to all lower case in-place.
 */
STR 
STR_ToLower(   STR      str );

/*! FUNCTION:  STR_StartsWith()
 *  SYNOPSIS:  Checks if start of <str> and <prefix> are equal up to the length of <prefix>.
 */
bool 
STR_StartsWith(   STR      str,
                  STR      prefix );

/*! FUNCTION:  STR_EndsWith()
 *  SYNOPSIS:  Checks if end of <str> and <suffix> are equal up to the length of <suffix>.
 */
bool 
STR_EndsWith(  STR      str,
               STR      suffix );

#endif /* _STR_EXT_H */
