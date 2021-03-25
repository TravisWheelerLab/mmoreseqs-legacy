/*******************************************************************************
 *  @file submat.h
 *  @brief SUBMAT Object (Substitution Matrix)
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

#ifndef _SUBMAT_H
#define _SUBMAT_H

/* Constructor */
SCORE_MATRIX* SCORE_MATRIX_Create();

/* Destructor */
void SCORE_MATRIX_Destroy( SCORE_MATRIX*   submat );

/* Construct SCORE_MATRIX object by parsing .submat file */
SCORE_MATRIX* SCORE_MATRIX_Load( char*    filename );

/* Set alphabet and the initialize score matrix based on size */
void SCORE_MATRIX_SetAlphabet( SCORE_MATRIX*   submat,
                                char*           alph );

/* Maps 2D-coords to 1D-coords in SUBSTITUTION MATRIX */
int SCORE_MATRIX_Keymap( SCORE_MATRIX*    submat,
                         char             q_ch, 
                         char             t_ch );

/* Get score from SCORE_MATRIX, given query/target chars. Returns reference.  */
float* SCORE_MATRIX_Score( SCORE_MATRIX*  submat, 
                           char           q_ch, 
                           char           t_ch );

/* Get score from SCORE_MATRIX, given query/target chars  */
float SCORE_MATRIX_GetScore( SCORE_MATRIX*  submat, 
                              char           q_ch, 
                              char           t_ch );

/* Output SCORE_MATRIX to FILE pointer */
void SCORE_MATRIX_Dump( SCORE_MATRIX*  submat, 
                        FILE*          fp );

#endif /* _SUBMAT_H */