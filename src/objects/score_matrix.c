/*******************************************************************************
 *  @file submat.c
 *  @brief SCORE_MATRIX Object.  
 *         Substitution Matrix for BLOSUM62, etc.
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/* local imports */
#include "../structs.h"
#include "../../utility.h"

/* header */
#include "score_matrix.h"

/* Constructor */
SCORE_MATRIX* SCORE_MATRIX_Create()
{
   SCORE_MATRIX*  submat;
   
   submat = (SCORE_MATRIX*) malloc( sizeof(SCORE_MATRIX) );
   if (submat == NULL) {
      perror("Error while malloc'ing SCORE_MATRIX.\n");
      exit(EXIT_FAILURE);
   }

   submat->scores = (float*) malloc( sizeof(float) * SUBMAT_SIZE );
   if (submat == NULL) {
      perror("Error while malloc'ing SCORES in SCORE_MATRIX.\n");
      exit(EXIT_FAILURE);
   }

   return submat;
}

/* Destructor */
void SCORE_MATRIX_Destroy(SCORE_MATRIX* submat)
{
   free(submat->filename);
   free(submat->scores);
   free(submat);
}

/* Parse .submat file and build SCORE_MATRIX object */
SCORE_MATRIX* SCORE_MATRIX_Load(char* _filename_)
{
   SCORE_MATRIX *submat;
   submat = SCORE_MATRIX_Create();
   submat->filename = strdup(_filename_);

   /* line reader objects */
   char*    line_buf       = NULL;
   size_t   line_buf_size  = 0;
   int      line_count     = 0;
   ssize_t  line_size      = 0;

   /* line parser objects */
   char     delim[]        = "\t";
   char*    parser         = NULL;
   int      x              = 0;
   int      y              = 0;
   int      key            = 0;
   float    score          = 0;
   char     a              = 0;
   char     b              = 0;

   /* open file */
   FILE *fp;
   fp = fopen(_filename_, "r");

   if (fp == NULL) {
      fprintf(stderr, "Error while opening file: %s\n", _filename_);
      exit(EXIT_FAILURE);
   }

   /* read a line */
   while ( ( line_size = getline ( &line_buf, &line_buf_size, fp ) ), line_size != -1 )
   {
      /* get next target character */
      b = AA[y];
      /* reset x at start of each line */
      x = 0;

      if (line_buf[0] != '#') /* if not a comment, parse line */
      {
         bool isNull = (parser != NULL);
         parser = strtok(line_buf, delim);

         // next row value is the first element on line
         if (parser != NULL) {
            b = parser[0];
            parser = strtok(NULL, delim);
         }
         while (parser != NULL) {
            /* get next query character */ 
            a = AA2[x]; 
            score = atof(parser);

            // check if cast score is valid
            if (!(score == 0 && parser[0] != '0')) {
               // map protein pair to int (both directions)
               *( SCORE_MATRIX_Score(submat,a,b) ) = score;
               *( SCORE_MATRIX_Score(submat,b,a) ) = score;
               x++;
            }
            // split on next tab
            parser = strtok(NULL, delim);
         }
      }
      y++;
   }
   fclose(fp);

   return submat;
}

/* Maps 2D-coords to 1D-coords in SUBSTITUTION MATRIX */
int SCORE_MATRIX_Keymap( char query_char, 
                         char target_char )
{
   int X = query_char - 'A';
   int Y = target_char - 'A';

   return (X * ALPHA_MAX) + Y;
}

/* Get score from SCORE_MATRIX, given query/target chars. Returns reference.  */
float* SCORE_MATRIX_Score( SCORE_MATRIX*  submat, 
                           char           q_ch, 
                           char           t_ch )
{
   int key = SCORE_MATRIX_Keymap(q_ch, t_ch);
   return &(submat->scores[key]);
}

/* Get score from SCORE_MATRIX, given query/target chars  */
float SCORE_MATRIX_Get_Score( SCORE_MATRIX*  submat, 
                              char           q_ch, 
                              char           t_ch )
{
   int key = SCORE_MATRIX_Keymap(q_ch, t_ch);
   return submat->scores[key];
}

/* Output SCORE_MATRIX to FILE pointer */
void SCORE_MATRIX_Dump( SCORE_MATRIX*  submat, 
                        FILE*          fp )
{
   fprintf(fp, "# ====== SCORE_MATRIX ======\n");
   fprintf(fp, "# FILENAME:\t%s\n", submat->filename);

   char query_char,target_char;
   float score;
   int i,j = 0;

   fprintf(fp, "#\t");
   for (i = 0; i < NUM_AMINO; i++)
   {
      fprintf(fp, "%c\t", AA[i]);
   }
   fprintf(fp, "\n");

   for (i = 0; i < NUM_AMINO; i++)
   {
      query_char = AA[i];
      fprintf(fp, "%c\t", query_char);
      for (j = 0; j < NUM_AMINO; j++)
      {
         target_char = AA[j];

         score = *( SCORE_MATRIX_Score(submat, query_char, target_char) );
         fprintf(fp, "%.3f\t", score);
      }
      fprintf(fp, "\n");
   }
   fprintf(fp, "# =========================\n\n");
}