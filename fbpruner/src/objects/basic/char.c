/*******************************************************************************
 *  FILE:      char.c
 *  PURPOSE:   CHAR Object
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
#include "../../utilities/_utilities.h"
#include "../_objects.h"

/* header */
#include "char.h"

/*
 *  FUNCTION:  CHAR_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
inline
char* CHAR_To_String( 	const CHAR   	d,
                     	char*       	buf )
{
   sprintf( buf, "%d", d );
   return buf;
}

/*
 *  FUNCTION:  CHAR_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
inline
int CHAR_Compare(  	const CHAR   a, 
                  	const CHAR   b )
{
   return (a - b);
}
