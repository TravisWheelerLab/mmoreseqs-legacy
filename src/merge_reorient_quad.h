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

int edgebounds_Merge_Reorient_Naive(EDGEBOUNDS *edg_fwd,
                                    EDGEBOUNDS *edg_bck,
                                    EDGEBOUNDS *edg_diag,
                                    EDGEBOUNDS *edg_row,
                                    int Q, int T,
                                    float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                                    float sp_MX[ NUM_NORMAL_STATES * (Q + 1) ]);

void edgebounds_Build_From_Cloud( EDGEBOUNDS*edg,
                                 int Q, int T,
                                 float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                                 int mode);

#endif /* _MERGE_REORIENT_QUAD_H */