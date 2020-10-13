/*******************************************************************************
 *  FILE:      str.c
 *  PURPOSE:   STR Object
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
#include "../../utilities/utilities.h"
#include "../objects.h"

/* header */
#include "str.h"

/*
 *  FUNCTION:  STR_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
inline
char* 
STR_To_String( const STR   d,
               char*       buf )
{
   sprintf( buf, "%s", d );
   return buf;
}

/*
 *  FUNCTION:  FLT_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
inline
int 
STR_Compare(  const STR   a, 
              const STR   b )
{
   return strcmp( a, b );
}
