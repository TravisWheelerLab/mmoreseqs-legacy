/*******************************************************************************
 *  FILE:      cloud_search_linear.h
 *  PURPOSE:   Cloud Search for Forward-Backward Pruning Algorithm. (LINEAR SPACE)
 *
 *  AUTHOR:    Dave Rich
 *  BUG:      
 *******************************************************************************/


#ifndef _CLOUD_SEARCH_LINEAR_ROWS_H
#define _CLOUD_SEARCH_LINEAR_ROWS_H

/* === FUNCTIONS === */

/* */
void cloud_Forward_Linear_Rows(  const SEQUENCE*      query, 
                                 const HMM_PROFILE*   target,
                                 const int            Q, 
                                 const int            T, 
                                 float*               st_MX, 
                                 float*               st_MX3,
                                 float*               sp_MX, 
                                 const ALIGNMENT*     tr,
                                 EDGEBOUNDS*          edg,
                                 const float          alpha, 
                                 const int            beta,
                                 const bool           test );

/* */
void cloud_Backward_Linear_Rows( const SEQUENCE*      query, 
                                 const HMM_PROFILE*   target,
                                 const int            Q, 
                                 const int            T, 
                                 float*               st_MX, 
                                 float*               st_MX3,
                                 float*               sp_MX, 
                                 const ALIGNMENT*     tr,
                                 EDGEBOUNDS*          edg,
                                 const float          alpha, 
                                 const int            beta,
                                 const bool           test );

#endif /* _CLOUD_SEARCH_LINEAR_ROWS_H */