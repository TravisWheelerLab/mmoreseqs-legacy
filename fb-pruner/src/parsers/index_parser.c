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
#include "structs.h"
#include "utilities.h"
#include "objects.h"

/* header */
#include "parsers.h"


/* ****************************************************************************************** *
 *  FUNCTION:  F_INDEX_Hmm_Build()
 *  SYNOPSIS:  Build F_INDEX object from .hmm file with list of hmm name and offset into file.
 *             
 *  METHOD:    [1]   scans .hmm file line-by-line
 *             [2]   if line is the header of a hmm model (line starts with "HMMER")
 *             [2a]  name of hmm and it offset location into a file added to index
 *
 *  ARGS:      <_filename_>   path to file to be indexed
 *
 *  RETURN:    F_INDEX object containing index of file
/* ****************************************************************************************** */
F_INDEX* F_INDEX_Hmm_Build( F_INDEX*    f_index,
                        const char* _filename_ )
{
   FILE*          fp             = NULL;
   F_INDEX_NODE   node;

   int            line_count     = 0;
   size_t         line_buf_size  = 0;
   ssize_t        line_size      = 0;
   char*          line_buf       = NULL;
   char*          token          = NULL;
   char*          delim          = "\t\n";

   long           prv_offset     = 0;
   long           cur_offset     = 0;

   int            id             = 0;
   char*          name           = NULL;
   int            mmseqs_id      = -1;

   /* create index object if necessary */
   if ( f_index == NULL ) {
      f_index = F_INDEX_Create();
   } else {
      F_INDEX_Reuse( f_index );
   }

   ERROR_free( f_index->source_path );
   f_index->source_path = strdup( _filename_ );
   /* index does not use mmseqs names by default */
   f_index->mmseqs_names = false;

   fp = fopen(_filename_, "r");
   if (fp == NULL) {
      fprintf(stderr, "ERROR: Unable to Open File '%s'\n", _filename_);
      exit(EXIT_FAILURE);
   }

   /* read file line-by-line */
   while( (line_size = getline(&line_buf, &line_buf_size, fp)), line_size >= 0 )
   {
      /* ignore commented lines */
      if (line_buf[0] == '#') continue;

      /* if starts new header, add to index */
      if ( strncmp(line_buf, "HMMER", 5)  == 0 ) 
      {
         cur_offset = prv_offset;

         while ( (line_size = getline(&line_buf, &line_buf_size, fp)), line_size >= 0 )
         {
            if ( strncmp(line_buf, "NAME", 4) == 0 ) 
            {
               int i = 0;
               for ( i = 4; line_buf[i] != ' '; i++) {}  /* skip whitespace after NAME */
               name = &line_buf[i];
               name[strlen(name)-1] = '\0';

               /* add to index */
               node.id        = id;
               node.name      = name;
               node.offset    = cur_offset;
               node.mmseqs_id = mmseqs_id;
               F_INDEX_Pushback( f_index, &node );

               id++;
               break;
            }
         }
      }
      
      line_count++;
      prv_offset = ftell(fp);
   }

   fclose(fp);
   ERROR_free(line_buf);
   return f_index;
}


/* ****************************************************************************************** *
 *  FUNCTION:  F_INDEX_Fasta_Build()
 *  SYNOPSIS:  Build F_INDEX object from .fasta file by fasta name and offset into file.
 *             
 *  METHOD:    [1]   scans .fasta file line-by-line
 *             [2]   if line is the header of a fasta sequence:
 *             [2a]  name of fasta and it offset location into a file added to index
 *
 *  ARGS:      <_filename_>   path to file to be indexed.
 *
 *  RETURN:    F_INDEX object containing index of file  
/* ****************************************************************************************** */
F_INDEX* F_INDEX_Fasta_Build( F_INDEX*    f_index,
                              const char* _filename_ )
{
   FILE*          fp             = NULL;
   F_INDEX_NODE   node;

   int            line_count     = 0;
   size_t         line_buf_size  = 0;
   ssize_t        line_size      = 0;
   char*          line_buf       = NULL;
   char*          header         = NULL;
   char*          token          = NULL;
   char*          delim          = "\t\n";

   long           prv_offset     = 0;
   long           cur_offset     = 0;

   int            id             = 0;
   char*          name           = NULL;
   int            mmseqs_id      = -1;

   /* create index object if necessary */
   if ( f_index == NULL ) {
      f_index = F_INDEX_Create();
   } else {
      F_INDEX_Reuse( f_index );
   }

   /* use filename as source path */
   ERROR_free( f_index->source_path );
   f_index->source_path = strdup( _filename_ );
   /* index does not use mmseqs names by default */
   f_index->mmseqs_names = false;

   /* open file */
   fp = fopen(_filename_, "r");
   if (fp == NULL) {
      fprintf(stderr, "ERROR: Unable to Open File => %s\n", _filename_);
      exit(EXIT_FAILURE);
   }

   /* read file line-by-line */
   while( (line_size = getline(&line_buf, &line_buf_size, fp)), line_size >= 0 )
   {
      /* ignore commented lines */
      if (line_buf[0] == '#') continue;

      /* remove end line */
      if(line_buf[line_size-1] == '\n') {
         line_buf[line_size-1] = '\0';
         line_size--;
      }

      /* if starts new header, add to index */
      if (line_buf[0] == '>') 
      {
         /* parse header */
         header = line_buf + 1;    /* skip '>' character */
         /* split header on spaces, get first element */
         token = strtok(header, " ");
         /* if name has structure: >db|id| */
         if ( strstr( token, "|" ) == NULL ) {
            name = strtok(token, "|");    /* 1st field = database */
            name = strtok(NULL, "|");    /* 2nd field = name */ 
         } 
         /* otherwise, just use the entire header */
         else {
            name = token;
         }

         node.id        = id;
         node.name      = name;
         node.offset    = prv_offset;
         node.mmseqs_id = -1;
         F_INDEX_Pushback( f_index, &node );
         id++;
      }
      
      line_count++;
      prv_offset = ftell(fp);
   }

   fclose(fp);
   ERROR_free(line_buf);
   return f_index;
}

/* ****************************************************************************************** *
 *  FUNCTION:  F_INDEX_Load()
 *  SYNOPSIS:  Build F_INDEX object from .idx file.
 *             
 *  METHOD:    [1]   scans .fasta file line-by-line
 *             [2]   line: {result_id} {offset} {}
 *
 *  ARGS:      <_filename_>   path to file to be indexed.
 *
 *  RETURN:    F_INDEX object containing index of file  
/* ****************************************************************************************** */
F_INDEX* F_INDEX_Load( F_INDEX*     f_index,
                       const char*  _filename_ )
{
   FILE*          fp             = NULL;
   F_INDEX_NODE   node;

   int            line_count     = 0;
   size_t         line_buf_size  = 0;
   ssize_t        line_size      = 0;
   char*          line_buf       = NULL;
   char*          token          = NULL;
   char*          delim          = "\t\n";

   long           prv_offset     = 0;
   long           cur_offset     = 0;

   int            id             = 0;
   char*          name           = NULL;
   int            mmseqs_id      = -1;

   /* create index object if necessary */
   if ( f_index == NULL ) {
      f_index = F_INDEX_Create();
   } else {
      F_INDEX_Reuse( f_index );
   }

   f_index->index_path = strdup( _filename_ );

   /* file open */
   fp = fopen(_filename_, "r");
   if (fp == NULL) {
      fprintf(stderr, "ERROR: Unable to Open File => %s\n", _filename_);
      exit(EXIT_FAILURE);
   }

   // /* parse header details */
   // while( (line_size = getline(&line_buf, &line_buf_size, fp)), line_size >= 0 )
   // {
   //    /* ignore commented lines */
   //    if (line_buf[0] == '#') continue;
   //    /* carot marks last line of header */
   //    if (line_buf[0] == '>') break;

   //    /* remove end line */
   //    if(line_buf[line_size-1] == '\n') {
   //       line_buf[line_size-1] = '\0';
   //       line_size--;
   //    }

   //    /* get next token, if empty line, skip it */
   //    token = strtok( line_buf, delim );
   //    if (token == NULL) continue;

   //    /* path to source */
   //    if ( strcmp( token, "SOURCE_PATH:" ) == 0 ) {
   //       token = strtok( NULL, delim );
   //       if ( strcmp( token, "(null)" ) != 0 ) {
   //          f_index->source_path = strdup( token );
   //       }

   //       line_count++;
   //       continue;
   //    }

   //    /* path to lookup */
   //    if ( strcmp( token, "LOOKUP_PATH:" ) == 0 ) {
   //       token = strtok( NULL, delim );
   //       if ( strcmp( token, "(null)" ) != 0 ) {
   //          f_index->lookup_path = strdup( token );
   //       }

   //       line_count++;
   //       continue;
   //    }

   //    /* number of nodes */
   //    if ( strcmp( token, "NUMBER_SEQS:" ) == 0 ) {
   //       token = strtok( NULL, delim );
   //       int size = atoi( token );
   //       F_INDEX_Resize( f_index, size + 1 );

   //       line_count++;
   //       continue;
   //    }

   //    /* whether using mmseqs names */
   //    if ( strcmp( token, "MMSEQS_NAMES:" ) == 0 ) {
   //       token = strtok( NULL, delim );
   //       if ( strcmp( token, "true" ) == 0 ) {
   //          f_index->mmseqs_names = true;
   //       } 
   //       else if ( strcmp( token, "false" ) == 0 ) {
   //          f_index->mmseqs_names = false;
   //       }

   //       line_count++;
   //       continue;
   //    }

   //    line_count++;
   // }

   /* read file line-by-line */
   while( (line_size = getline(&line_buf, &line_buf_size, fp)), line_size >= 0 )
   {
      /* ignore commented lines */
      if (line_buf[0] == '#') continue;

      /* remove end line */
      if(line_buf[line_size-1] == '\n') {
         line_buf[line_size-1] = '\0';
         line_size--;
      }

      /* first token is the id */
      token = strtok( line_buf, delim );
      if (token == NULL) continue;
      id = atoi( token );

      /* second token is the offset */
      if ( (token = strtok( NULL, delim )), token == NULL) continue;
      cur_offset = atol( token ); 

      /* third token is the name */
      if ( (token = strtok( NULL, delim )), token == NULL) continue;
      name = token;

      node.id        = id;
      node.name      = name;
      node.offset    = cur_offset;
      node.mmseqs_id = -1;
      F_INDEX_Pushback( f_index, &node );

      line_count++;
   }

   fclose(fp);
   ERROR_free(line_buf);
   return f_index;
}

/* update f_index using mmseqs lookup table */
void F_INDEX_Lookup_Update( F_INDEX*   f_index, 
                            char*      _lookup_filepath_ ) 
{
   FILE*          fp             = NULL;
   F_INDEX_NODE*  node;

   int            line_count     = 0;
   size_t         line_buf_size  = 0;
   ssize_t        line_size      = 0;
   char*          line_buf       = NULL;
   char*          token          = NULL;
   char*          delim          = "\t\n";

   long           prv_offset     = 0;
   long           cur_offset     = 0;

   int            id             = 0;
   char*          name           = NULL;
   int            mmseqs_id      = -1;

   /* add meta data */
   f_index->lookup_path = strdup( _lookup_filepath_ );

   /* open file */
   fp = fopen( _lookup_filepath_, "r" );
   if (fp == NULL) {
      fprintf(stderr, "ERROR: Unable to Open File => %s\n", _lookup_filepath_ );
      exit(EXIT_FAILURE);
   }

   /* parse header details */
   while( (line_size = getline(&line_buf, &line_buf_size, fp)), line_size >= 0 )
   {
      /* ignore commented lines */
      if (line_buf[0] == '#') continue;

      /* first token is id */
      token = strtok(line_buf, delim);
      id = atoi(token);

      /* get the idth node */
      node = &(f_index->nodes[id]);

      /* second token is name */
      token = strtok(line_buf, delim);
      name = token;

      /* update name to match mmseqs lookup */
      ERROR_free(node->name);
      node->name = strdup( name );
   }

   fclose(fp);
   ERROR_free(line_buf);
}

/* load index file */
F_INDEX* F_INDEX_Plus_Load(   F_INDEX*       f_index,
                              const char*    _filename_ )
{
   FILE*          fp             = NULL;
   F_INDEX_NODE   node;

   int            line_count     = 0;
   size_t         line_buf_size  = 0;
   ssize_t        line_size      = 0;
   char*          line_buf       = NULL;
   char*          token          = NULL;
   char*          delim          = "\t\n";

   long           prv_offset     = 0;
   long           cur_offset     = 0;

   int            id             = 0;
   char*          name           = NULL;
   int            mmseqs_id      = -1;

   /* create index object if necessary */
   if ( f_index == NULL ) {
      f_index = F_INDEX_Create();
   } else {
      F_INDEX_Reuse( f_index );
   }

   f_index->index_path = strdup( _filename_ );

   /* file open */
   fp = fopen(_filename_, "r");
   if (fp == NULL) {
      fprintf(stderr, "ERROR: Unable to Open File => %s\n", _filename_);
      exit(EXIT_FAILURE);
   }

   /* read file line-by-line */
   while( (line_size = getline(&line_buf, &line_buf_size, fp)), line_size >= 0 )
   {
      /* ignore commented lines */
      if (line_buf[0] == '#') continue;

      /* remove end line */
      if(line_buf[line_size-1] == '\n') {
         line_buf[line_size-1] = '\0';
         line_size--;
      }

      /* first token is the id */
      token = strtok( line_buf, delim );
      if (token == NULL) continue;
      id = atoi( token );

      /* second token is the offset */
      if ( (token = strtok( NULL, delim )), token == NULL) continue;
      cur_offset = atol( token ); 

      /* first token is the id */
      token = strtok( line_buf, delim );
      if (token == NULL) continue;
      mmseqs_id = atoi( token );

      /* third token is the name */
      if ( (token = strtok( NULL, delim )), token == NULL) continue;
      name = token;

      F_INDEX_Pushback( f_index, &node );

      line_count++;
   }

   fclose(fp);
   ERROR_free(line_buf);
   return f_index;
}