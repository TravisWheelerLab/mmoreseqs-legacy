/*******************************************************************************
 *  FILE:      _vectors.h
 *  PURPOSE:   Vectors of Primitive Datatypes
 *
 *  AUTHOR:    Dave Rich
 *   NOTES:
 *******************************************************************************/

#ifndef _VECTOR_H
#define _VECTOR_H

/* this is just used for generating the other basic datatype files */
#include "vector_template.h"

/* primitives (wrappers) */
#include "vector_char.h"
#include "vector_double.h"
#include "vector_float.h"
#include "vector_int.h"

/* composite (more than a single primitive) */
#include "vector_range.h"
#include "vector_trace.h"
#include "vector_bound.h"
// #include "vectors_str.h"

/* generic datatype (can potentially be dependent on any basic type) */
// #include "vector_gen.h"

#endif /* _VECTOR_H */
