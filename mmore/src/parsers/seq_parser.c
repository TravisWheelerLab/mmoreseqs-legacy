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
#include "../objects/structs.h"
#include "../utilities/_utilities.h"
#include "../objects/_objects.h"

/* header */
#include "_parsers.h"

/* parse .fasta file and build SEQUENCE object */
void SEQUENCE_Fasta_Parse( SEQUENCE*   seq,
                           char*       filename,
                           long        offset )
{
   /* parser vars */
   READER*  reader         = NULL;
   STR      line           = NULL;
   int      line_count     = -1;       /* line counter of current line in file */
   char*    line_buf       = NULL;     /* pointer to start of buffered line */
   char*    header         = NULL;     /* header pointer */
   char*    name           = NULL;     /* name pointer */
   char*    token          = NULL;     /* token for splitting string */
   size_t   line_buf_size  = 0;        /* length of entire <line_buf> array */
   size_t   line_size      = 0;        /* length of current line in <line_buf> array */
   int      lines_read     = 0;
   int      num_seqs       = 0;
   int      max_seqs       = 1;        
   int      seq_len        = 0;

   /* clear pre-existing data */
   SEQUENCE_Reuse( seq );
   SEQUENCE_SetTextfield( &seq->filename, filename );

   /* open file reader */
   reader = READER_Create( filename );
   READER_Open( reader );
   READER_JumpTo( reader, offset );

   /* read file line-by-line */
   while ( READER_NextLine( reader ), READER_Is_EndOfFile( reader ) == false )
   {
      line_count++;
      bool is_eof    = READER_Is_EndOfFile( reader );
      line           = READER_GetLine( reader );
      
      /* check if line is a header */
      if (line[0] == '>')
      {
         /* if onto next seq, return (only captures first sequence) */
         if ( num_seqs > 0 ) {
            break;
         }
         /* increment number of sequences captured */
         num_seqs += 1;

         /* omit first ">" character */
         header = BUFFER_Move( reader->buffer, 1 );

         /* first capture whole header line */
         SEQUENCE_SetTextfield( &(seq->header), header );

         /* split header on spaces, get first field delimited by " " */
         token = strtok(header, " ");
         /* if name has structure: >db|id|, then get the third field delimited by "|" */
         if ( strstr( token, "|" ) != NULL ) {
            name = strtok(token, "|");
            name = strtok(NULL, "|");
         } 
         /* otherwise, just use the whole field */
         else {
            name = token;
         }
         SEQUENCE_SetTextfield( &(seq->name), name ); 
      }
      else /* otherwise, append to current sequence */
      {
         /* capitalize all characters (model doesn't account for sequence uncertainty) */
         line = STR_ToUpper( line );
         /* append to growing sequence */
         SEQUENCE_Append_Seq( seq, line );
      }
   }

   READER_Close( reader );
   READER_Destroy( reader );
}
