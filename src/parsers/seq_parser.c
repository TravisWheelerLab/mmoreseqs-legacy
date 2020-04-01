/*******************************************************************************
 *  FILE:      seq_parser.h
 *  PURPOSE:   Parses .fasta files into SEQUENCE object
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
#include <math.h>
#include <ctype.h>

/* local imports */
#include "objects/structs.h"
#include "objects/hmm_profile.h"
#include "objects/sequence.h"
#include "utility.h"

/* header */
#include "hmm_parser.h"

/* parse .fasta file and build SEQUENCE object */
void SEQUENCE_Fasta_Parse( SEQUENCE*   seq,
                           char*       _filename_,
                           long        offset )
{
   /* parser vars */
   FILE*    fp             = NULL;    
   int      line_count     = -1;       /* line counter of current line in file */
   char*    line_buf       = NULL;     /* pointer to start of buffered line */
   size_t   line_buf_size  = 0;        /* length of entire <line_buf> array */
   size_t   line_size      = 0;        /* length of current line in <line_buf> array */

   int      num_seqs       = 0;        
   int      seq_len        = 0;

   SEQUENCE_Reuse(seq);
   SEQUENCE_Set_Textfield(&seq->filename, _filename_);

   /* open file */
   fp = fopen(_filename_, "r");
   /* check for file read error */
   if (fp == NULL)
   {
      char*  str = NULL;
      fprintf(stderr, "ERROR: Bad FILE POINTER for SEQUENCE PARSER => %s\n", _filename_ );
      exit(EXIT_FAILURE);
   }
   fseek(fp, offset, SEEK_SET);

   /* read file line-by-line */
   while ( ( line_size = getline ( &line_buf, &line_buf_size, fp ) ), line_size != -1 )
   {
      line_count++;

      /* remove newline from end of line */
      if( line_buf[line_size-1] == '\n' ) {
         line_buf[line_size-1] = '\0';
         line_size -= 1;
      }

      /* check if line is a header */
      if (line_buf[0] == '>')
      {
         num_seqs++;
         SEQUENCE_Set_Textfield(&seq->name, line_buf);

         /* if onto next seq, reset counters */
         /* NOTE: Currently only accepts one SEQUENCE */
         if (num_seqs > 1) {
            break;
         }

         char* line_ptr = line_buf;
         char* token    = NULL;
         int   i        = 0;
         
         while( token = strtok_r(line_ptr, "/|\n", &line_ptr) )
         {
            // printf("[%d] %s\n", i, token);
            i++;
         }
      }
      else /* otherwise, append to current sequence */
      {
         /* capitalize all characters */
         for (int i=0; i<line_size; i++) {
            line_buf[i] = toupper(line_buf[i]);
         }

         SEQUENCE_Append_Seq(seq, line_buf);
      }
   }

   /* free line buffer */
   free(line_buf);
   /* close file */
   fclose(fp);
}
