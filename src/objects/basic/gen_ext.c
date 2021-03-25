/*******************************************************************************
 *  FILE:      gen_ext.h
 *  PURPOSE:   GEN Object.  A union which can hold most primitive datatypes.
 *             Extended functions.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *    - None.
 *  TODO:
 *    - Could create an faster (though unsafe) ToString() function.
 *    - Add support for int RANGE type.
 *    - Need to implement Compare().
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
#include "_basic.h"
#include "gen.h"

/* lookup for datatype sizes */
const size_t DATATYPE_SIZES[] = {
   0,                /* DATATYPE_NONE */
   sizeof(int),      /* DATATYPE_INT */
   sizeof(float),    /* DATATYPE_FLOAT */
   sizeof(float),    /* DATATYPE_FLOAT_EXP */
   sizeof(double),   /* DATATYPE_DOUBLE */
   sizeof(double),   /* DATATYPE_DOUBLE_EXP */
   sizeof(long),     /* DATATYPE_LONG */
   sizeof(char*),    /* DATATYPE_STRING */
   sizeof(char),     /* DATATYPE_CHAR */
   sizeof(bool),     /* DATATYPE_BOOL */
   sizeof(void*)     /* DATATYPE_POINTER */
};

/*! FUNCTION:  GEN_Wrap()
 *  SYNOPSIS:  Wrap into a GEN struct. 
 *             <data> should be reference to the data to be stored.
 *             <type> should be one of the enumerated datatypes.
 *             <size> should be sizeof(<type>). 
 *    RETURN:  GEN struct.
 */
inline
GEN 
GEN_Wrap( 	const void*       data,
            const DATATYPE    type,
            const size_t      size )
{
   GEN      gdata;
   size_t   gsize;

   gdata.type = type;
   /* should resist compiler trying to cast <data> */
   memcpy( &gdata.data, data, size );

   return gdata;
}

/*! FUNCTION:  GEN_GetSize()
 *  SYNOPSIS:  Get the size of an enumerated datatype.
 *    RETURN:  Pointer to <buf>, or NULL if error.
 */
inline
size_t 
GEN_GetSize( 	const DATATYPE    type )
{
   return DATATYPE_SIZES[type];
}

