/*******************************************************************************
 *  FILE:      objects.h
 *  PURPOSE:   All /object/ folder headers.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *    - None Known.
 *  NOTES:
 *    - This include folder needs to respect inter-object folder dependencies.
 *    - The plan is to extract most struct definitions and place them in the files with their function defs.
 *******************************************************************************/

#ifndef _OBJECTS_H
#define _OBJECTS_H

/* declares all datatypes */
#include "structs.h"

/* basic types (no dependencies) */
#include "basic/_basic.h"

/* vector types (only dependent on basic types) */
#include "vectors/_vectors.h"

/* normal matrix types */
#include "matrix/_matrix.h"

/* sparse matrix */
#include "matrix_sparse/_matrix_sparse.h"

/* complex objects (can be dependent on basics, vectors, and matrices) */
#include "alignment.h"
#include "arg_opts.h"
#include "args.h"
#include "clock.h"
#include "debugger.h"
#include "domain_def.h"
#include "f_index.h"
#include "hmm_profile.h"
#include "hmm_bg.h"
#include "mystring.h"
#include "results.h"
#include "m8_results.h"
#include "score_matrix.h"
#include "sequence.h"
#include "scriptrunner.h"
#include "commandline.h"

/* worker (dependent on most object types) */
#include "worker_thread.h"
#include "worker.h"

/* io (dependent on all so that it can read in / write any data type) */
#include "io/_io.h"

#endif /* _OBJECTS_H */