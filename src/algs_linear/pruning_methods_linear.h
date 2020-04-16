/*******************************************************************************
 *  FILE:      pruning_methods.h
 *  PURPOSE:   Assorted Pruning Methods for Cloud Search
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

/* ****************************************************************************************** *
 *  
 *  FUNCTION:  Prune_Xdrop_Trim_Edges.
 *  SYNOPSIS:  
 *
 *  ARGS:      <edg>          List of edgebounds (stored by anti-diags)
 *             <d_cnt>        Current antidiagonal being pruned
 *             <alpha>        Pruning ratio
 *             <beta>         Number of edgebounds to traverse before pruning
 *             <st_MX3>       DP Matrix
 *             <d0>           The current anti-diagonal (mod-mapped to linear matrix)
 *             <d1>           The preveding anti-diagonal (mod-mapped to linear matrix)
 *             <lb_1>         The preceding anti-diagonal's 
 *
 *  RETURN:    No return, but <edg> is updated with latest diag edgebounds.
 *
/* ****************************************************************************************** */
void Prune_Xdrop_Trim3( EDGEBOUNDS*  edg,
                        int          beta,
                        float        alpha,
                        MATRIX_3D*   st_MX3,
                        int          d_cnt,
                        int          d0,
                        int          d1 );