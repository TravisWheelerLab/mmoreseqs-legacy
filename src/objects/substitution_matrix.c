/*******************************************************************************
 *  @file submat.c
 *  @brief SUBSTITUTION_MATRIX Object
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
#include "substitution_matrix.h"

/* Constructor */
SUBSTITUTION_MATRIX* SUBSTITUTION_MATRIX_Create()
{
   SUBSTITUTION_MATRIX *submat;
   submat = malloc( sizeof(SUBSTITUTION_MATRIX) );
   if (submat == NULL) 
   {
      perror("Error while malloc'ing SUBSTITUTION_MATRIX.\n");
      exit(EXIT_FAILURE);
   }
   submat->scores = (float *)malloc( sizeof(float) * SUBMAT_SIZE );
   if (submat == NULL) 
   {
      perror("Error while malloc'ing SCORES in SUBSTITUTION_MATRIX.\n");
      exit(EXIT_FAILURE);
   }

   return submat;
}

/* Destructor */
void SUBSTITUTION_MATRIX_Destroy(SUBSTITUTION_MATRIX *submat)
{
   free(submat->filename);
   free(submat->scores);
   free(submat);
}

/* Parse .submat file and build SUBSTITUTION_MATRIX object */
SUBSTITUTION_MATRIX* SUBSTITUTION_MATRIX_Load(char *_filename_)
{
   SUBSTITUTION_MATRIX *submat;
   submat = SUBSTITUTION_MATRIX_Create();

   /* line reader objects */
   char        *line_buf = NULL;
   size_t      line_buf_size = 0;
   int         line_count = 0;
   ssize_t     line_size;

   /* line parser objects */
   char        delim[] = "\t";
   char        *parser;
   int         x,y;
   int         key;
   float       score;
   char        a,b;
   x = 0; y = 0;

   /* open file */
   FILE *fp;
   fp = fopen(_filename_, "r");

   if (fp == NULL) 
   {
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
            a = AA[x]; 
            score = atof(parser);

            // check if cast score is valid
            if (!(score == 0 && parser[0] != '0')) {
               // map protein pair to int (both directions)
               key = SUBSTITUTION_MATRIX_Keymap(a,b);
               submat->scores[key] = score;
               key = SUBSTITUTION_MATRIX_Keymap(b,a);
               submat->scores[key] = score;
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
int SUBSTITUTION_MATRIX_Keymap(char query_char, 
                               char target_char)
{
   int X = query_char - 'A';
   int Y = target_char - 'A';

   return (X * ALPHA_MAX) + Y;
}

/* Get score from SUBSTITUTION_MATRIX */
float SUBSTITUTION_MATRIX_Get_Score(SUBSTITUTION_MATRIX *submat, 
                                    char query_char, 
                                    char target_char)
{
   int key = SUBSTITUTION_MATRIX_Keymap(query_char, target_char);
   return submat->scores[key];
}

/* Output SUBSTITUTION_MATRIX to FILE pointer */
void SUBSTITUTION_MATRIX_Dump(SUBSTITUTION_MATRIX *submat, 
                              FILE *fp)
{
   fprintf(fp, "\n");
   fprintf(fp, "====== SUBSTITUTION_MATRIX ======\n");

   char query_char,target_char;
   float score;
   int i,j = 0;

   fprintf(fp, "\t");
   for (i = 0; i < NUM_AMINO; i++)
   {
      fprintf(fp, "[%c]\t", AA[i]);
   }
   fprintf(fp, "\n");

   for (i = 0; i < NUM_AMINO; i++)
   {
      query_char = AA[i];
      fprintf(fp, "[%c]\t", query_char);
      for (j = 0; j < NUM_AMINO; j++)
      {
         target_char = AA[j];

         score = SUBSTITUTION_MATRIX_Get_Score(submat, query_char, target_char);
         fprintf(fp, "%.3f\t", score);
      }
      fprintf(fp, "\n");
   }
   fprintf(fp, "====================\n\n");
}