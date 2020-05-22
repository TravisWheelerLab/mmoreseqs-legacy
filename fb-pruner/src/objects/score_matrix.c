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
#include "structs.h"
#include "utilities.h"
#include "objects.h"

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

   submat->filename  = NULL;
   submat->alph      = NULL;
   submat->scores    = NULL;

   SCORE_MATRIX_Set_Alphabet( submat, ALPH_AMINO_CHARS );

   return submat;
}

/* Destructor */
void SCORE_MATRIX_Destroy(SCORE_MATRIX* submat)
{
   free(submat->filename);
   free(submat->alph);
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
   char*    token          = NULL;

   /* line parser objects */
   bool     hasHeader      = false;
   char*    header         = strdup( ALPH_AMINO_CHARS );
   char     delim[]        = " \t\n";
   char*    parser         = NULL;

   int      x              = 0;  /* index of row in matrix */
   int      y              = 0;  /* index of column in matrix */
   char     a              = 0;  /* char at row index */
   char     b              = 0;  /* char at column index */
   float    score          = 0;  /* match score */

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
      /* get first non-whitespace word */
      parser = strtok(line_buf, delim);

      /* if comment or blank line, skip */
      if ( parser == NULL || parser[0] == '#' ) {
         continue;
      }

      /* header leads with a > character */
      if ( parser[0] == '>' && !hasHeader )
      {
         header = (char*) realloc( header, sizeof(char) * 256 );
         int i = 0;

         /* parse header elements */
         while ( (parser = strtok(NULL, delim) ), parser != NULL )  
         {
            header[i] = parser[0];
            i++;
         }
         header[i] = '\0';
         hasHeader = true;
         SCORE_MATRIX_Set_Alphabet( submat, header );
         continue;
      }

      /* if line begins with letter from alphabet, we are on a data row */
      if ( strchr( header, parser[0] ) != NULL ) 
      {
         /* first element on line is member of alphabet */
         b = parser[0];

         /* tab-delimited scores */
         x = 0;
         while ( (parser = strtok(NULL, delim) ), parser != NULL ) 
         {
            /* get next query character */ 
            a = header[x]; 
            score = atof(parser);

            // map protein pair to int (both directions)
            *( SCORE_MATRIX_Score( submat, a, b ) ) = score;
            *( SCORE_MATRIX_Score( submat, b, a ) ) = score;

            x++;
         }
         y++;
         continue;
      }

   }
   fclose(fp);

   free(header);
   free(line_buf);

   return submat;
}

/* Set alphabet and the initialize score matrix based on size */
void SCORE_MATRIX_Set_Alphabet( SCORE_MATRIX*   submat,
                                char*           alph )
{
   int alph_len      = strlen( alph );
   submat->alph      = NULL;
   submat->alph      = strdup( alph );
   submat->alph_len  = alph_len;
   
   /* initialize map */
   for ( int i = 0; i < 256; i++ ) {
      submat->map[i] = -1;
   }

   /* create map from char -> index */
   for ( int i = 0; i < strlen(alph); i++ ) {
      submat->map[alph[i]] = i;
   }

   submat->scores = (float*) realloc( submat->scores, sizeof(float) * alph_len * alph_len );
   if (submat == NULL) {
      perror("Error while malloc'ing SCORES in SCORE_MATRIX.\n");
      exit(EXIT_FAILURE);
   }
}

/* Maps 2D-coords to 1D-coords in SUBSTITUTION MATRIX */
int SCORE_MATRIX_Keymap( SCORE_MATRIX*    submat,
                         char             q_ch, 
                         char             t_ch )
{
   int X = submat->map[q_ch];
   int Y = submat->map[t_ch];

   return ( (X * submat->alph_len) + Y );
}

/* Get score from SCORE_MATRIX, given query/target chars. Returns reference.  */
inline
float* SCORE_MATRIX_Score( SCORE_MATRIX*  submat, 
                           char           q_ch, 
                           char           t_ch )
{
   int key = SCORE_MATRIX_Keymap(submat, q_ch, t_ch);
   return &(submat->scores[key]);
}

/* Get score from SCORE_MATRIX, given query/target chars  */
inline
float SCORE_MATRIX_Get_Score( SCORE_MATRIX*  submat, 
                              char           q_ch, 
                              char           t_ch )
{
   int key = SCORE_MATRIX_Keymap(submat, q_ch, t_ch);
   return submat->scores[key];
}

/* Output SCORE_MATRIX to FILE pointer */
void SCORE_MATRIX_Dump( SCORE_MATRIX*  submat, 
                        FILE*          fp )
{
   fprintf(fp, "# ====== SCORE_MATRIX ======\n");
   fprintf(fp, "# FILE:\t%s\n", submat->filename);
   fprintf(fp, "# ALPH:\t%s\n", submat->alph);

   char     q_ch     = 0;
   char     t_ch     = 0;
   float    score    = 0.0;

   int      i        = 0; 
   int      j        = 0;
   char*    alph     = submat->alph;
   int      alph_len = submat->alph_len;

   /* print header */
   fprintf(fp, ">\t");
   for (i = 0; i < alph_len; i++) 
   {
      fprintf(fp, "%c\t", alph[i]);
   }
   fprintf(fp, "\n");

   /* print each row */
   for (i = 0; i < alph_len; i++) 
   {
      fprintf(fp, "%c\t", alph[i]);
      for (j = 0; j < alph_len; j++)
      {
         score = *( SCORE_MATRIX_Score( submat, alph[i], alph[j] ) );
         fprintf(fp, "%f\t", score);
      }
      fprintf(fp, "\n");
   }
   fprintf(fp, "# ============================\n");
}

/* for now, just create generic amino alphabet */
ALPHABET* ALPHABET_Create( ALPHABET alph_type )
{

}