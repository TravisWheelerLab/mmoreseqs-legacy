/*******************************************************************************
 *  FILE:      reader.c
 *  PURPOSE:   READER Class. Helps with reading from files.
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
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

/* local imports */
#include "structs.h"
#include "../utilities/utilities.h"
#include "objects.h"

/* header */
#include "debugger.h"

/*  FUNCTION:  READER_Create()
 *  SYNOPSIS:
 */
READER* READER_Create( char* filename )
{
	READER* reader       = NULL;
   reader               = ERROR_malloc( sizeof(READER) );
   reader->filename     = strdup( filename );
   /* buffer data */
   reader->buffer       = NULL;
   reader->N            = 0;
   reader->Nalloc       = 0;
   reader->token        = NULL;
   reader->tokens       = VECTOR_CHAR_Create();
   /* location in file */
   reader->cur_offset   = 0;
   reader->file_size    = 0;
   reader->lines_read   = 0;
   reader->is_eof       = false;

   READER_Open( reader );

   return reader;
}

/*  FUNCTION:  READER_Destroy()
 *  SYNOPSIS:
 */
void* READER_Destroy( READER* reader )
{
   if ( reader == NULL ) return NULL;

   READER_Close( reader );

   VECTOR_CHAR_Destroy( reader->tokens );
   ERROR_free( reader->filename );
   ERROR_free( reader->buffer );
   ERROR_free( reader );

   return NULL;
}

/*  FUNCTION:  READER_Open()
 *  SYNOPSIS:
 */
int READER_Open( READER* reader )
{
   reader->fp = ERROR_fopen( reader->filename, "r" );

   return STATUS_SUCCESS;
}

/*  FUNCTION:  READER_Close()
 *  SYNOPSIS:
 */
int READER_Close( READER* reader )
{
   ERROR_fclose( reader->fp );
   reader->fp = NULL;

   return STATUS_SUCCESS;
}

/*  FUNCTION:  READER_Rewind()
 *  SYNOPSIS:  Set file pointer to beginning of file.
 */
int READER_Rewind( READER* reader )
{
   int status = fseek( reader->fp, 0, SEEK_SET );
   return status;
}

/*  FUNCTION:  READER_JumpTo()
 *  SYNOPSIS:  Jump to nth byte of into file pointer.
 */
int READER_JumpTo(   READER*     reader, 
                     long int    offset )
{
   int status = fseek( reader->fp, offset, SEEK_SET );
   return status;
}

/*  FUNCTION:  READER_GetLine()
 *  SYNOPSIS:
 */
char* READER_GetLine( READER* reader )
{
   if ( reader->is_eof ) return NULL;

   reader->N = getline( &reader->buffer, &reader->Nalloc, reader->fp );

    /* if at normal line */
   if ( reader->N != -1 ) {
      return reader->buffer;
   } 
   /* if at end of file */
   else {
      reader->is_eof = true;
      return NULL;
   }
}

/*  FUNCTION:  READER_Split()
 *  SYNOPSIS:  Split line currently in buffer into tokens.
 */
int READER_SplitLine(   READER*  reader, 
                        char*    delimiter )
{
   VECTOR_CHAR_Reuse( reader->tokens );

   return reader->tokens;
}