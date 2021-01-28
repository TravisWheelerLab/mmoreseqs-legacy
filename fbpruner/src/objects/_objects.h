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

/* io */
#include "io/reader.h"

/* basic datatypes / wrappers */
#include "basic/_basic.h"

/* vectors */
#include "vectors/vector_bound.h"
#include "vectors/vector_char.h"
#include "vectors/vector_double.h"
#include "vectors/vector_float.h"
#include "vectors/vector_int.h"
#include "vectors/vector_range.h"
#include "vectors/vector_trace.h"
#include "vectors/vector_template.h"

/* matrix */
#include "matrix/matrix_2d.h"
#include "matrix/matrix_3d.h"
#include "matrix/matrix_3d_sparse.h"

/* objects */
#include "alignment.h"
#include "args.h"
#include "clock.h"
#include "debugger.h"
#include "domain_def.h"
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

#endif /* _OBJECTS_H */