/*******************************************************************************
 *  FILE:      _basic.h
 *  PURPOSE:   Basic data functionality.
 *  AUTHOR:    Dave Rich
 * 
 *  DESC:      API for building vectors, matrices, and maps.
 *             For simple non-dynamic types, most of these are pass-through functions that should be optimized out by compiler.
 *             Promises the following API for each datatype (XXX is a template datatype):
 *                XXX      XXX_Create( XXX data ) 
 *                XXX      XXX_Destroy( XXX data )
 *                XXX      XXX_Set( XXX data )
 *                size_t   XXX_Sizeof( XXX data )
 *                int      XXX_Compare( XXX a, XXX b )
 *                bool     XXX_Equals( XXX a, XXX b )
 *                void     XXX_Swap( XXX* a, XXX* b )
 *                STR      XXX_To_String( XXX data, STR buffer ) 
 *                -  NOTE: Need to add a size parameter to _To_String(), so buffer doesn't overflows.
 *
 *             Files with "xxx_ext.h" stores functions outside core API (may be moved to objects/ later).
 *
 *   NOTES:
 *    - There is one level of potential dependency: composite types can incorporate primitive types.
 *******************************************************************************/

#ifndef _BASIC_H
#define _BASIC_H

/* special case: because _To_String() requires knowledge of STR, that goes first */
#include "str.h"

/* this is just used for generating the other basic datatype files */
#include "template.h"

/* primitives (wrappers) */
#include "bool.h"
#include "char.h"
#include "double.h"
#include "float.h"
#include "int.h"

/* composite datatypes (more than a single primitive) */
#include "bound.h"
#include "range.h"
#include "trace.h"

/* dynamically allocated / pointers datatypes */
#include "str.h"
#include "str_ext.h"
#include "xstr.h"
#include "ptr.h"

/* generic datatype (can potentially be dependent on any basic type) */
#include "gen.h"
#include "gen_ext.h"

#endif /* _BASIC_H */
