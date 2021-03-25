/*******************************************************************************
 *  FILE:      xstr.h
 *  PURPOSE:   XSTR object.
 *             String builder class. Better for manipulating strings than STR.
 *
 *  AUTHOR:    Dave Rich
 *  NOTES:
 *    - WIP.       
 *******************************************************************************/

#ifndef _X_STRING_H
#define _X_STRING_H

/*! FUNCTION:  X_STRING_Create()
 *  SYNOPSIS:  Allocates data to hold <xstr> with length of <chars>
 */
 XSTR* 
 XSTR_Create( char* chars );

 /*! FUNCTION:  XSTR_Destroy()
 *  SYNOPSIS:  Destroys <str> data and frees memory.
 */
XSTR* 
XSTR_Destroy( XSTR* str );

/*! FUNCTION:  XSTR_Copy()
 *  SYNOPSIS:  Creates a deep copy of <src>, stores it in <dest> and returns it.
 */
XSTR* 
XSTR_Copy(  XSTR*    src,
            XSTR*    dest );

/*! FUNCTION:  XSTR_SetSize()
 *  SYNOPSIS:  Allocates <str> to store a string of length <L>.
 *             
 */
XSTR* 
XSTR_SetSize( XSTR*      str,
               size_t     L );

/*! FUNCTION:  X_STRING_Set()
 *  SYNOPSIS:  Set <str> to hold <src>.
 */
XSTR* 
XSTR_Set( XSTR*   str,
          char*   src );

/*! FUNCTION:  XSTR_Concat()
 *  SYNOPSIS:  Copies <str_1> and <str_2> (in order) into new <str>.
 *             Returns <str>.
 */
XSTR* 
XSTR_Concat( XSTR*   str_1,
             XSTR*   str_2 );

/*! FUNCTION:  XSTR_Append()
 *  SYNOPSIS:  Appends <str_add> <str> and <str_2> (in order) into new <str>.
 *             Returns <str>.
 */
XSTR* 
XSTR_Append( XSTR*   str,
             XSTR*   str_add );

#endif /* _X_STRING_H */