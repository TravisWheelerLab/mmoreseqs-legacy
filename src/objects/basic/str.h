/*******************************************************************************
 *  FILE:      str.c
 *  PURPOSE:   STR Object ( wraps *char )
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _STR_H
#define _STR_H

/*! FUNCTION:  STR_Create()
 *  SYNOPSIS:  Create <str>, fills with input <str_chrs> and return pointer to <str>.
 */
STR
STR_Create(  const char*   chrs );

/*! FUNCTION:  STR_Destroy()
 *  SYNOPSIS:  Destroys <str>, frees memory and returns NULL pointer.
 */
STR
STR_Destroy( STR   str );

/*! FUNCTION:  STR_Clear()
 *  SYNOPSIS:  Clear <data>.  If pointer data, sets to null. Otherwise, do nothing. 
 */
STR
STR_Clear( STR   data );

/*! FUNCTION:  STR_ToString()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *    RETURN:  Pointer to <buf>
 */
char* 
STR_ToString(  const STR   data,
               char*       buf );

/*! FUNCTION:  STR_FromString()
 *  SYNOPSIS:  Extracts data from string.
 *
 *    RETURN:  Pointer to <buf>
 */
STR
STR_FromString(  char*   str );

/*! FUNCTION:  STR_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *    RETURN:  POS if (a > b), 
 *             0 if equal, 
 *             NEG if (a < b)
 */
int 
STR_Compare(  const STR   a, 
              const STR   b );

/*! FUNCTION:  STR_CompareTo()
 *  SYNOPSIS:  Generic compare. Casts then compares <a> and <b>.
 *    RETURN:  POS if (a > b), 
 *             0 if equal, 
 *             NEG if (a < b)
 */
int 
STR_CompareTo(    const void*   a, 
                  const void*   b );

#endif /* _STR_H */
