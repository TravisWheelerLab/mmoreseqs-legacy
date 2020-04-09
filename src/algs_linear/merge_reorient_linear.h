/*******************************************************************************
 *  FILE:      merge_reorient_linear.c
 *  PURPOSE:   Functions for merging multiple EDGEBOUND objects and 
 *             reorienting from diagonal-wise to row-wise.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/


#ifndef _MERGE_REORIENT_LINEAR_H
#define _MERGE_REORIENT_LINEAR_H

/* === INCLUDES === */
// #include "objects/structs.c"
// #include "objects/edgebound.c"

/* === FUNCTIONS === */
void edgebounds_Reflect(EDGEBOUNDS *edg);

EDGEBOUNDS* EDGEBOUNDS_Merge(const int         Q, 
                             const int         T,
                             const EDGEBOUNDS* edg_1,
                             const EDGEBOUNDS* edg_2);

EDGEBOUNDS* EDGEBOUNDS_Reorient(const int         Q, 
                                const int         T,
                                const EDGEBOUNDS* edg_in);

#endif /* _MERGE_REORIENT_LINEAR_H */