/*******************************************************************************
 *  @file merge_reorient_linear.c
 *  @brief  Functions for merging multiple EDGEBOUND objects and reorienting from antidiagonal-wise to row-wise.
 *
 *  @synopsis
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/


#ifndef _MERGE_REORIENT_LINEAR_H
#define _MERGE_REORIENT_LINEAR_H

void edgebounds_Reflect(EDGEBOUNDS *edg);

void edgebounds_Merge(int Q, int T,
                     EDGEBOUNDS *edg_fwd,
                     EDGEBOUNDS *edg_bck,
                     EDGEBOUNDS *edg_new);

void edgebounds_Reorient(int Q, int T,
                         EDGEBOUNDS *edg_diag,
                         EDGEBOUNDS *edg_row);

#endif /* _MERGE_REORIENT_LINEAR_H */