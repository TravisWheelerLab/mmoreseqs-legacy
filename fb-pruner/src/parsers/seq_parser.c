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
#include "structs.h"
#include "utilities.h"
#include "objects.h"

/* header */
#include "parsers.h"

/* parse .fasta file and build SEQUENCE object */
void SEQUENCE_Fasta_Parse( SEQUENCE*   seq,
                           char*       _filename_,
                           long        offset )
{
   /* parser vars */
   FILE*    fp             = NULL;    
   int      line_count     = -1;       /* line counter of current line in file */
   char*    line_buf       = NULL;     /* pointer to start of buffered line */
   char*    header         = NULL;     /* header pointer */
   char*    name           = NULL;     /* name pointer */
   char*    token          = NULL;     /* token for splitting string */
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
         /* if onto next seq, return (only captures first sequence) */
         if (num_seqs > 0) {
            break;
         }

         num_seqs++;
         /* omit ">" */
         header = line_buf + 1;

         /* first capture whole header line */
         SEQUENCE_Set_Textfield(&seq->header, header);

         /* split header on spaces, get first element */
         token = strtok(header, " ");
         /* if name has structure: >db|id| */
         if ( strstr( token, "|" ) == NULL ) {
            name = strtok(token, "|");
            name = strtok(NULL, "|");
         } 
         /* otherwise, just use the  */
         else {
            name = token;
         }
         SEQUENCE_Set_Textfield(&seq->name, name); 
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
   ERROR_free(line_buf);
   /* close file */
   fclose(fp);
}
