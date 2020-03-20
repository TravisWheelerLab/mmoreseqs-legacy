/*******************************************************************************
 *  FILE:      index_parser.c
 *  PURPOSE:   Build an INDEX object from: .idx, .hmm, or .fasta file
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

/* local imports */
#include "../structs.h"
#include "objects/f_index.h"
#include "objects/vectors/vector_int.h"

/* header */
#include "index_parser.h"


/* ****************************************************************************************** *
 *  
 *  FUNCTION:  F_INDEX_Hmm_Build
 *  SYNOPSIS:  Build F_INDEX object from .hmm file with list of hmm name and offset into file.
 *             
 *  METHOD:    [1]   scans .hmm file line-by-line
 *             [2]   if line is the header of a hmm model (line starts with "HMMER")
 *             [2a]  name of hmm and it offset location into a file added to index
 *
 *  ARGS:      <_filename_>   path to file to be indexed
 *
 *  RETURN:    F_INDEX object containing index of file
 *
/* ****************************************************************************************** */
F_INDEX* F_INDEX_Hmm_Build( const char*    _filename_ )
{
   F_INDEX*    f_index        = NULL;
   FILE*       fp             = NULL;

   char*       line_buf       = NULL;
   size_t      line_buf_size  = 0;

   ssize_t     line_size      = 0;
   int         line_count     = 0;

   long        prv_offset     = 0;
   long        cur_offset     = 0;
   char*       name           = NULL;

   /* first pass */
   printf("FILE: %s\n", _filename_);
   f_index = F_INDEX_Create(_filename_);
   fp      = fopen(_filename_, "r");

   if (fp == NULL) {
      fprintf(stderr, "ERROR: Unable to Open File => %s\n", _filename_);
      exit(EXIT_FAILURE);
   }

   /* read file line-by-line */
   while( (line_size = getline(&line_buf, &line_buf_size, fp)), line_size >= 0 )
   {
      /* if starts new header, add to index */
      if ( strncmp(line_buf, "HMMER", 5) ) 
      {
         cur_offset = prv_offset;

         while ( (line_size = getline(&line_buf, &line_buf_size, fp)), line_size >= 0 )
         {
            if ( strncmp(line_buf, "NAME", 4) ) {
               int i = 0;
               for ( i = 4; line_buf[i] != ' '; i++) {}
               name = &line_buf[i];
            }
            F_INDEX_PushBack( f_index, (F_INDEX_NODE){line_buf, cur_offset} );
         }
      }
      
      line_count++;
      prv_offset = ftell(fp);
      // printf("%s [%d=>%lu]:{%d/%d} %s\n", _filename_, line_count, tell, line_size, line_buf_size, line_buf);
   }

   free(line_buf);
   return f_index;
}


/* ****************************************************************************************** *
 *  
 *  FUNCTION:  F_INDEX_Fasta_Build
 *  SYNOPSIS:  Build F_INDEX object from .fasta file by fasta name and offset into file.
 *             
 *  METHOD:    [1]   scans .fasta file line-by-line
 *             [2]   if line is the header of a fasta sequence:
 *             [2a]  name of fasta and it offset location into a file added to index
 *
 *  ARGS:      <_filename_>   path to file to be indexed.
 *
 *  RETURN:    F_INDEX object containing index of file  
 *
/* ****************************************************************************************** */
F_INDEX* F_INDEX_Fasta_Build( const char*    _filename_ )
{
   F_INDEX*    f_index        = NULL;
   FILE*       fp             = NULL;

   char*       line_buf       = NULL;
   size_t      line_buf_size  = 0;

   ssize_t     line_size      = 0;
   int         line_count     = 0;

   long        cur_offset     = 0;
   long        prv_offset     = 0;
   char*       name           = NULL;

   /* first pass */
   printf("FILE: %s\n", _filename_);
   f_index = F_INDEX_Create(_filename_);
   fp      = fopen(_filename_, "r");

   if (fp == NULL) {
      fprintf(stderr, "ERROR: Unable to Open File => %s\n", _filename_);
      exit(EXIT_FAILURE);
   }

   /* read file line-by-line */
   while( (line_size = getline(&line_buf, &line_buf_size, fp)), line_size >= 0 )
   {
      /* remove end line */
      if(line_buf[line_size-1] == '\n') {
         line_buf[line_size-1] = '\0';
         line_size--;
      }

      /* if starts new header, add to index */
      if (line_buf[0] == '>') {
         F_INDEX_PushBack( f_index, (F_INDEX_NODE){line_buf, prv_offset} );
      }
      
      line_count++;
      prv_offset = ftell(fp);
      // printf("%s [%d=>%lu]:{%d/%d} %s\n", _filename_, line_count, tell, line_size, line_buf_size, line_buf);
   }

   free(line_buf);
   return f_index;
}

/* load index file */
F_INDEX* F_INDEX_Load( const char*   _filename_ )
{
   F_INDEX*    f_index        = NULL;
   FILE*       fp             = NULL;

   char*       line_buf       = NULL;
   size_t      line_buf_size  = 0;

   ssize_t     line_size      = 0;
   int         line_count     = 0;

   long        cur_offset     = 0;
   long        prv_offset     = 0;
   char*       name           = NULL;

   /* first pass */
   printf("FILE: %s\n", _filename_);
   f_index = F_INDEX_Create(_filename_);
   fp      = fopen(_filename_, "r");

   if (fp == NULL) {
      fprintf(stderr, "ERROR: Unable to Open File => %s\n", _filename_);
      exit(EXIT_FAILURE);
   }

   /* parse header details */


   /* read file line-by-line */
   while( (line_size = getline(&line_buf, &line_buf_size, fp)), line_size >= 0 )
   {
      /* remove end line */
      if(line_buf[line_size-1] == '\n') {
         line_buf[line_size-1] = '\0';
         line_size--;
      }

      /* */


      line_count++;
   }

   free(line_buf);
   return f_index;
}