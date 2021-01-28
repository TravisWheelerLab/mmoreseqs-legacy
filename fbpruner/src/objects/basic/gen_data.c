/*******************************************************************************
 *  FILE:      gen_data.c
 *  PURPOSE:   GEN_DATA Object. A union that can store most primitive data types.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *    - None.
 *  TODO:
 *    - Could create an faster (though unsafe) To_String() function
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
#include "bound.h"

/*! FUNCTION:  GEN_DATA_Create()
 *  SYNOPSIS:  Create a GEN_DATA struct. 
 *             <data> should be reference to the data to be stored.
 *             <type> should be one of the enumerated datatypes.
 *             <size> should be sizeof(<type>). 
 *    RETURN:  Pointer to GEN_DATA struct.
 */
inline
GEN_DATA 
GEN_DATA_Create( 	const void*       data,
                  const DATATYPE    type,
                  const size_t      size )
{
   GEN_DATA gdata;
   gdata.type = type;
   /* should resist compiler trying to cast <data> */
   memcpy( &gdata.data, data, size );

   return gdata;
}

/*! FUNCTION:  GEN_DATA_To_String()
 *  SYNOPSIS:  Create a string representation of <data>.
 *             If it is of float-like type, formats with <sig_digits> as number of significant digits.
 *             Stores it in a char* buffer <buf>. Caller must have preallocated buffer. 
 *             Strings are truncated to buf_size.
 *    RETURN:  Pointer to <buf>, or NULL if error.
 */
inline
char* 
GEN_DATA_To_String( 	const GEN_DATA    data,
                     char*             buf,
                     const int         buf_size,
                     const int         sig_digits )
{
   /* data array */
   if ( data.type == DATATYPE_NONE ) {
      snprintf( buf, buf_size, "%s\t", "--" );
   }
   elif ( data.type == DATATYPE_STRING ) {
      snprintf( buf, buf_size, "%s\t", data.data.s );
   }
   elif ( data.type == DATATYPE_FLOAT ) {
      snprintf( buf, buf_size, "%.*f\t", sig_digits, data.data.f );
   }
   elif ( data.type == DATATYPE_BOOL ) {
      snprintf( buf, buf_size, "%s\t", (data.data.b == true ? "T" : "F") );
   }
   elif ( data.type == DATATYPE_INT ) {
      snprintf( buf, buf_size, "%d\t", data.data.i );
   }
   elif ( data.type == DATATYPE_LONG ) {
      snprintf( buf, buf_size, "%ld\t", data.data.l );
   }
   elif ( data.type == DATATYPE_CHAR ) {
      snprintf( buf, buf_size, "%c\t", data.data.c );
   }
   elif ( data.type == DATATYPE_FLOAT_EXP ) {
      snprintf( buf, buf_size, "%.*e\t", sig_digits, data.data.f );
   }
   elif ( data.type == DATATYPE_DOUBLE ) {
      snprintf( buf, buf_size, "%.*f\t", sig_digits, data.data.d );
   }
   elif ( data.type == DATATYPE_DOUBLE_EXP ) {
      snprintf( buf, buf_size, "%.*e\t", sig_digits, data.data.d );
   }
   elif ( data.type == DATATYPE_POINTER ) {
      snprintf( buf, buf_size, "%p\t", data.data.p );
   }
   else {
      fprintf( stderr, "ERROR: Invalid data type! VALID RANGE = (%d,%d), VALUE GIVEN = %d\n",
         1, NUM_DATATYPES, data.type );
      ERRORCHECK_exit(EXIT_FAILURE);
      return NULL;
   }
}

/*! FUNCTION:  GEN_DATA_Compare()
 *  SYNOPSIS:  Compare <a> and <b>.
 *             Sort first according to type, then by value.
 *    RETURN:  Positive if (a > b), 
 *             Zero if equal, 
 *             Negative if (a < b)
 */
inline
int 
GEN_DATA_Compare(    const GEN_DATA   a, 
                     const GEN_DATA   b )
{
   
}
