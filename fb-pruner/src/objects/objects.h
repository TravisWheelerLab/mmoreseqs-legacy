/*******************************************************************************
 *  FILE:      objects.h
 *  PURPOSE:   All /object/ folder headers.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _OBJECTS_H
#define _OBJECTS_H

/* declares all datatypes */
#include "structs.h"

/* basic datatypes / wrappers */
#include "bound.h"
#include "char.h"
#include "double.h"
#include "float.h"
#include "int.h"
#include "range.h"
#include "trace.h"
#include "template.h"

/* vectors */
#include "vector_bound.h"
#include "vector_char.h"
#include "vector_double.h"
#include "vector_float.h"
#include "vector_int.h"
#include "vector_range.h"
#include "vector_trace.h"
#include "vector_template.h"

/* matrix */
#include "matrix_2d.h"
#include "matrix_3d.h"
#include "matrix_3d_sparse.h"

/* objects */
#include "alignment.h"
#include "args.h"
#include "clock.h"
#include "debugger.h"
#include "edgebound.h"
#include "edgebound_rows.h"
#include "f_index.h"
#include "hmm_profile.h"
#include "hmm_bg.h"
#include "mystring.h"
#include "results.h"
#include "score_matrix.h"
#include "sequence.h"
#include "worker.h"
#include "worker_thread.h"
#include "x_string.h"
#include "reader.h"

#endif /* _OBJECTS_H */