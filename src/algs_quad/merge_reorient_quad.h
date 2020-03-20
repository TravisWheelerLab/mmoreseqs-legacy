/*******************************************************************************
 *  @file merge_reorient_quad.h
 *  @brief  Functions for merging multiple EDGEBOUND objects and reorienting from antidiagonal-wise to row-wise.
 *
 *  @synopsis
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/


#ifndef _MERGE_REORIENT_QUAD_H
#define _MERGE_REORIENT_QUAD_H

/* === INCLUDES === */
// #include "objects/structs.h"
// #include "objects/edgebound.h"

/* === FUNCTIONS === */
int EDGEBOUNDS_Merge_Reorient_Naive(const EDGEBOUNDS* edg_fwd,
                                    const EDGEBOUNDS* edg_bck,
                                    EDGEBOUNDS*       edg_diag,
                                    EDGEBOUNDS*       edg_row,
                                    const int         Q, 
                                    const int         T,
                                    float*            st_MX,
                                    float*            sp_MX);

void EDGEBOUNDS_Build_From_Cloud(EDGEBOUNDS* edg,
                                 const int   Q, 
                                 const int   T,
                                 float*      st_MX,
                                 int         mode);

#endif /* _MERGE_REORIENT_QUAD_H */