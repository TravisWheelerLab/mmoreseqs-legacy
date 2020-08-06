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

/* matrix */
#include "matrix_2d.h"
#include "matrix_3d.h"
#include "matrix_3d_sparse.h"

/* vectors */
#include "vector_bound.h"
#include "vector_float.h"
#include "vector_int.h"
#include "vector_range.h"
#include "vector_range_2d.h"
#include "vector_trace.h"
#include "vector_char.h"

/* objects */
#include "alignment.h"
#include "args.h"
#include "clock.h"
#include "edgebound.h"
#include "edgebound_rows.h"
#include "f_index.h"
#include "hmm_profile.h"
#include "mystring.h"
#include "results.h"
#include "score_matrix.h"
#include "sequence.h"
#include "worker.h"

#endif /* _RESULTS_H */