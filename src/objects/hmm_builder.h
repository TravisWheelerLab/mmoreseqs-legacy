/*******************************************************************************
 *  @file hmm.c
 *  @brief HMM Object
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

#ifndef _HMM_BUILDER_H
#define _HMM_BUILDER_H

/* constructor */
HMM* HMM_BUILDER_Create();

/* destructor */
void HMM_BUILDER_Destory( HMM*  hmm );

#endif /* _HMM_BUILDER_H */