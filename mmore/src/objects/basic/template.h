/*******************************************************************************
 *  FILE:      template.h
 *  PURPOSE:   XXX Object
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _XXX_H
#define _XXX_H

/* include where XXX is defined */
#include "../structs.h"

/*! FUNCTION:  XXX_Create()
 *  SYNOPSIS:  Create an instance of <data>.
 */
XXX
XXX_Create( const XXX   data );

/*! FUNCTION:  XXX_Destroy()
 *  SYNOPSIS:  Destroy instance of <data>.
 */
XXX
XXX_Destroy( XXX   data );

/*! FUNCTION:  XXX_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a preallocated char* buffer <buf> of length <buf_size>.
 *    RETURN:  Pointer to <buf>.
 */
char* 
XXX_To_String(    const XXX   data,
                  char*       buf );

/*! FUNCTION:  XXX_From_String()
 *  SYNOPSIS:  Interpret string representation in char* buffer <buf>.
 *             Stores result in <data>.
 *    RETURN:  <STATUS_SUCCESS> if no errors.
 */
int 
XXX_From_String(  const char*    buf,
                  XXX            data );

/*! FUNCTION:  XXX_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
int 
XXX_Compare(   const XXX   a, 
               const XXX   b );

/*! FUNCTION:  XXX_Swap()
 *  SYNOPSIS:  Swap values of <a> and <b>
 */
void
XXX_Swap(   XXX*    a,
            XXX*    b );

#endif /* _XXX_H */
