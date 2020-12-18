/*******************************************************************************
 *  FILE:      bound.c
 *  PURPOSE:   BOUND Object
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
#include "bound.h"

/*
 *  FUNCTION:  BOUND_To_String()
 *  SYNOPSIS:  Create a string representation of data <d>.
 *             Stores it in a char* buffer <buf>.
 *
 *    RETURN:  Pointer to <buf>
 */
inline
char* BOUND_To_String( 	const BOUND   d,
                     	char*       buf )
{
   sprintf( buf, "{ [%d] (%d,%d) }", d.id, d.lb, d.rb );
   return buf;
}

/*
 *  FUNCTION:  BOUND_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *
 *    RETURN:  pos if (a > b), 
 *             0 if equal, 
 *             neg if (a < b)
 */
inline
int BOUND_Compare(   const BOUND   a, 
                     const BOUND   b )
{
   if ( (a.id - b.id) != 0 ) {
      return (a.id - b.id);
   }

   if ( (a.lb - b.lb) != 0 ) {
      return (a.lb - b.lb);
   }

   return (a.rb - b.rb);
}
