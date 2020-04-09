/*******************************************************************************
 *  FILE:      cloud_search_linear.c
 *  PURPOSE:   Cloud Search for Forward-Backward Pruning Alg. (LINEAR SPACE)
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/


#ifndef _CLOUD_SEARCH_LINEAR_H
#define _CLOUD_SEARCH_LINEAR_H

/* ------------------------------------------------------------------------------------------ *
 *
 *  FUNCTION: cloud_Forward_Linear()
 *  SYNOPSIS: 
 *
 *  ARGS:      <query>        Query sequence, 
 *             <target>       HMM model,
 *             <Q>            Query length, 
 *             <T>            Target length,
 *             <st_MX3>       Normal State (Match, Insert, Delete) Matrix (Linear Space),
 *             <st_MX>        Normal State (Match, Insert, Delete) Matrix (Quadratic Space),
 *             <sp_MX>        Special State (J,N,B,C,E) Matrix
 *
 *  RETURN:    None.
 *
 * ------------------------------------------------------------------------------------------- */
void cloud_Forward_Linear( const SEQUENCE*    query, 
                           const HMM_PROFILE* target,
                           const int          Q, 
                           const int          T, 
                           MATRIX_3D*         st_MX, 
                           MATRIX_3D*         st_MX3,
                           MATRIX_2D*         sp_MX, 
                           const ALIGNMENT*   tr,
                           EDGEBOUNDS*        edg,
                           const float        alpha, 
                           const int          beta,
                           const bool         test);

/* ------------------------------------------------------------------------------------------ *
 *
 *  FUNCTION: cloud_Backward_Linear()
 *  SYNOPSIS: 
 *
 *  ARGS:      <query>        Query sequence, 
 *             <target>       HMM model,
 *             <Q>            Query length, 
 *             <T>            Target length,
 *             <st_MX3>       Normal State (Match, Insert, Delete) Matrix (Linear Space),
 *             <st_MX>        Normal State (Match, Insert, Delete) Matrix (Quadratic Space),
 *             <sp_MX>        Special State (J,N,B,C,E) Matrix
 *
 *  RETURN:    Returns the final score of the Forward Algorithm.
 *
 * ------------------------------------------------------------------------------------------- */
void cloud_Backward_Linear(const SEQUENCE*    query, 
                           const HMM_PROFILE* target,
                           const int          Q, 
                           const int          T, 
                           MATRIX_3D*         st_MX, 
                           MATRIX_3D*         st_MX3,
                           MATRIX_2D*         sp_MX, 
                           const ALIGNMENT*   tr,
                           EDGEBOUNDS*        edg,
                           const float        alpha, 
                           const int          beta,
                           const bool         test);

#endif /* _CLOUD_SEARCH_LINEAR_H */