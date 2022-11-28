/*******************************************************************************
 *  - FILE:      str.c
 *  - DESC:    STR Object ( wraps *char )
 *******************************************************************************/

#ifndef _STR_H
#define _STR_H

/*! FUNCTION:  STR_Create()
 *  SYNOPSIS:  Create <str>, fills with input <str_chrs> and return pointer to
 * <str>.
 */
STR STR_Create(const char* chrs);

/*! FUNCTION:  STR_Destroy()
 *  SYNOPSIS:  Destroys <str>, frees memory and returns NULL pointer.
 */
STR STR_Destroy(STR str);

/*! FUNCTION:  STR_Clear()
 *  SYNOPSIS:  Clear <data>.  If pointer data, sets to null. Otherwise, do
 * nothing.
 */
STR STR_Clear(STR data);

/*! FUNCTION:  STR_ToString()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *    RETURN:  Pointer to <buf>
 */
char* STR_ToString(const STR data, char* buf);

/*! FUNCTION:  STR_FromString()
 *  SYNOPSIS:  Extracts data from string.
 *
 *    RETURN:  Pointer to <buf>
 */
STR STR_FromString(char* str);



#endif /* _STR_H */
