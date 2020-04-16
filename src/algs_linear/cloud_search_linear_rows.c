/*******************************************************************************
 *  FILE:      cloud_search_linear.h
 *  PURPOSE:   Cloud Search for Forward-Backward Pruning Algorithm. (LINEAR SPACE)
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/* data stuctures and utility functions */
#include "objects/structs.h"
#include "utilities/utility.h"
#include "testing.h"

/* objects */
#include "objects/edgebound.h"
#include "objects/hmm_profile.h"
#include "objects/sequence.h"
#include "objects/alignment.h"
#include "objects/vectors/vector_range_2d.h"

/* file parsers */
#include "hmm_parser.h"

/* header */
#include "cloud_search_linear_rows.h"

/* 
 *       NOTE: CONVERSION - row-coords => diag-coords
 *       MX(i-1, j-1) => MX3(d_2, k-1)
 *       MX(i  , j-1) => MX3(d_1, k  ) 
 *       MX(i-1, j  ) => MX3(d_1, k-1)
 *
 *       MX(i+1, j+1) => MX3(d_2, k+1)
 *       MX(i  , j+1) => MX3(d_1, k  )
 *       MX(i+1, j  ) => MX3(d_1, k+1)
 */


