/*******************************************************************************
 *  FILE:      work_sparse_mx.h
 *  PURPOSE:   Workflow Subroutines for Posterior Algorithms.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _WORK_SPARSE_MX_H
#define _WORK_SPARSE_MX_H

/*! FUNCTION:  	WORK_build_sparse_matrix()
 *  SYNOPSIS:  	Builds <st_SMX_fwd> and <st_SMX_bck> based on the <edg_rows> shape. 
 */
void 
WORK_build_sparse_matrix( WORKER*  worker );


#endif /* _WORK_SPARSE_MX_H */