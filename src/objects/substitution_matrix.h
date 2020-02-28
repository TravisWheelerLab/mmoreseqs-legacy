/*******************************************************************************
 *  @file submat.h
 *  @brief SUBMAT Object (Substitution Matrix)
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

#ifndef _SUBMAT_H
#define _SUBMAT_H

/* === INCLUDES === */
#include <stdio.h>

/* === CONSTANTS === */
/* max ASCII value of x - 'A' for alphabet */
#define ALPHA_MAX 26
#define SUBMAT_SIZE ALPHA_MAX*ALPHA_MAX

/* === OBJECTS === */
typedef struct {
   char *filename;
   float *scores;
} SUBSTITUTION_MATRIX;

/* === FUNCTIONS === */
/* constructor */
SUBSTITUTION_MATRIX* SUBSTITUTION_MATRIX_Create();
/* Destructor */
void SUBSTITUTION_MATRIX_Destroy(SUBSTITUTION_MATRIX *submat);
/* Parse .submat file and build SUBSTITUTION_MATRIX object */
SUBSTITUTION_MATRIX* SUBSTITUTION_MATRIX_Load(char *_filename_);

/* Get SUBSTITUTION_MATRIX Score, given query and target character */
float SUBSTITUTION_MATRIX_Get_Score(SUBSTITUTION_MATRIX *submat, 
                                    char query_char, 
                                    char target_char);
/* Maps 2D-coords to 1D-coords in SUBSTITUTION MATRIX */
int SUBSTITUTION_MATRIX_Keymap(char query_char, 
                               char target_char);

/* Output SUBSTITUTION_MATRIX to FILE pointer */
void SUBSTITUTION_MATRIX_Dump(SUBSTITUTION_MATRIX *submat, 
                              FILE *fp);

#endif /* _SUBMAT_H */