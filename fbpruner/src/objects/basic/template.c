/*******************************************************************************
 *  FILE:       template.c
 *  PURPOSE:    XXX Object
 *
 *  AUTHOR:     Dave Rich
 * 	NOTE: 		Generated using TEMPLATE via "scripts/builder-helper/build_basic_wrapper_files.sh"
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
#include "template.h"

/*
 *  FUNCTION:  XXX_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
inline
char* XXX_To_String( const XXX   d,
                     char*       buf )
{
   sprintf( buf, "%d", d );
   return buf;
}

/*
 *  FUNCTION:  XXX_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
inline
int XXX_Compare(  const XXX   a, 
                  const XXX   b )
{
   return (a - b);
}