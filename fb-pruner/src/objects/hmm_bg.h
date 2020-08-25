/*******************************************************************************
 *  FILE:      hmm_bg.c
 *  PURPOSE:   HMM_BG Object.
 *
 *  AUTHOR:    Dave Rich
 *  BUGS:  
 *******************************************************************************/

#ifndef _HMM_BG_H
#define _HMM_BG_H

/*
 *  FUNCTION:  HMM_BG_Create()
 *  SYNOPSIS:  
 */
HMM_BG* HMM_BG_Create();

/*
 *  FUNCTION:  HMM_BG_Destroy()
 *  SYNOPSIS:
 */
void* HMM_BG_Destroy( HMM_BG* bg );

/* FUNCTION:  HMM_BG_SetSequence()
 * SYNOPSIS:  Set the sequence to create digitized sequence.
 */
void HMM_BG_SetSequence( 	HMM_BG*		bg, 
									SEQUENCE* 	seq );

/*
 *  FUNCTION:  HMM_BG_SetLength()
 *  SYNOPSIS:  
 */
void HMM_BG_SetLength( 	HMM_BG*		bg, 
								int 			L );

/*
 *  FUNCTION:  HMM_BG_SetLength()
 *  SYNOPSIS:  
 */
void HMM_BG_NullOne( 	const HMM_BG* 		bg,
								int 					L, 
								float*				null_sc );

/*
 *  FUNCTION:  HMM_BG_SetFilter()
 *  SYNOPSIS:  
 */
void HMM_BG_SetFilter(	HMM_BG* 			bg, 
								int 				M, 
								const float*	compo );

/*
 *  FUNCTION:  HMM_BG_FilterScore()
 *  SYNOPSIS:  
 */
void HMM_BG_FilterScore( 	HMM_BG* 			bg,
                           SEQUENCE* 		q_seq,
                           float*			result );



#endif /* _HMM_BG_H */
