/*******************************************************************************
 *  FILE:      _basic.h
 *  PURPOSE:   Basic data type wrappers (simple compare, to_string, from_string operators).
 *
 *  AUTHOR:    Dave Rich
 *   NOTES:
 *    - There is one level of potential dependency: composite types can incorporate primitive types.
 *******************************************************************************/

#ifndef _BASIC_H
#define _BASIC_H

/* this is just used for generating the other basic datatype files */
#include "template.h"

/* primitives (wrappers) */
#include "bool.h"
#include "char.h"
#include "double.h"
#include "float.h"
#include "int.h"

/* composite (more than a single primitive) */
#include "bound.h"
#include "range.h"
#include "trace.h"
#include "str.h"
#include "str_ext.h"
#include "xstr.h"

/* generic datatype (can potentially be dependent on any basic type) */
#include "gen.h"
#include "gen_ext.h"

#endif /* _BASIC_H */
