/*******************************************************************************
 *  - FILE:      _vectors.h
 *  - DESC:    Vectors of Basic Datatypes.
 *             Promises the following API for each datatype (XXX is a template
 *datatype): VEC_XXX*    VEC_XXX_Create( (?) )
 *                   - TODO: should make _Create() take a pointer arg so that it
 *is consistent with basic types. VEC_XXX*    VEC_XXX_Destroy( VEC_XXX vec )
 *                size_t      VEC_XXX_Sizeof( VEC_XXX vec )
 *                int         VEC_XXX_Compare( VEC_XXX a, VEC_XXX b )
 *                int         VEC_XXX_ToString( VEC_XXX, STR buffer )
 *
 *                VEC_XXX*    VEC_XXX_Create_by_Size( size_t size )
 *                VEC_XXX*    VEC_XXX_WrapArray( XXX* array, size_t size )
 *                VEC_XXX*    VEC_XXX_UnwrapArray( VEC_XXX* vec )
 *                size_t      VEC_XXX_GetSize( VEC_XXX vec )
 *                void        VEC_XXX_Resize( VEC_XXX vec, size_t size )
 *                void        VEC_XXX_GrowTo( VEC_XXX vec )
 *                void        VEC_XXX_Reuse( VEC_XXX vec )
 *                XXX         VEC_XXX_Get( VEC_XXX vec, size_t i )
 *                XXX*        VEC_XXX_GetX( VEC_XXX vec, size_t i )
 *                void        VEC_XXX_Set( VEC_XXX vec, size_t, XXX val )
 *                void        VEC_XXX_Swap( VEC_XXX vec, size_t i, size_t j )
 *                int         VEC_XXX_Push( VEC_XXX vec, XXX val )
 *                int         VEC_XXX_Pushback( VEC_XXX vec, XXX val )
 *                XXX         VEC_XXX_Pop( VEC_XXX vec )
 *                XXX         VEC_XXX_Popback( VEC_XXX vec )
 *                int         VEC_XXX_Sort( VEC_XXX vec )
 *                int         VEC_XXX_Search( VEC_XXX vec, XXX val )
 *                int         VEC_XXX_SearchFirst( VEC_XXX vec, XXX val )
 *                int         VEC_XXX_SearchLast( VEC_XXX vec, XXX val )
 *
 *******************************************************************************/

#ifndef _VECTOR_H
#define _VECTOR_H

/* this is just used for generating the other basic datatype files */
#include "vector_template.h"

/* primitive datatypes (aliased by simple typedef macro) */
#include "vector_char.h"
#include "vector_double.h"
#include "vector_float.h"
#include "vector_int.h"

/* composite datatypes (more than a single primitive) */
#include "vector_bound.h"
#include "vector_range.h"
#include "vector_trace.h"

/* dynamically allocated / pointers datatypes */
#include "vector_ptr.h"
#include "vector_str.h"

/* XSTR is the string builder class - functionally it is a VEC_CHAR */
// #include "vector_xstr.h"

/* generic datatype (can potentially be dependent on any basic type) */
// #include "vector_gen.h"

#endif /* _VECTOR_H */
