/*******************************************************************************
 *  FILE:      system_io.c
 *  PURPOSE:   Tools for file input/output and other OS-related functions.
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
#include <limits.h>

/* local imports */
#include "structs.h"
#include "_objects.h"

/* header */
#include "_utilities.h"
#include "system_io.h"

/*! FUNCTION:  	SYSTEMIO_Exists()
 *  SYNOPSIS:  	Checks whether <filename> exists.
 */
bool
SYSTEMIO_Exists( const STR filename )
{
   FILE* file;
   if ( file = fopen( filename, "r" ) ) {
      fclose( file );
      return true;
   }
   else {
      return false;
   }
}

/*! FUNCTION:  	SYSTEMIO_HasWritePermission()
 *  SYNOPSIS:  	Checks whether user has permission to write to <filename>.
 */
bool
SYSTEMIO_HasWritePermission( const STR filename )
{
   if ( access( filename, W_OK ) == 0 ) {
      return true;
   }
   else {
      return false;
   }
}

/*! FUNCTION:  	SYSTEMIO_HasReadPermission()
 *  SYNOPSIS:  	Returns whether user has permission to read from <filename>.
 */
bool
SYSTEMIO_HasReadPermission( const STR filename )
{
   if ( access( filename, R_OK ) == 0 ) {
      return true;
   }
   else {
      return false;
   }
}

/*! FUNCTION:  	SYSTEMIO_GetDirectory()
 *  SYNOPSIS:  	Returns current working directory.
 *                Malloc's string data if passed NULL pointer.
 */
STR 
SYSTEMIO_GetDirectory( STR old_str )
{
   STR   new_str;
   char  cwd[PATH_MAX + 1];
   
   /* clean up old data */
   old_str = STR_Destroy( old_str );

   /* try to get current working directory */
   if ( getcwd(cwd, sizeof(cwd)) != NULL ) {
      new_str = STR_Create( cwd );
      return new_str;
   } else {
      fprintf( stderr, "ERROR: Could not get current working directory.\n");
      return NULL;
   }
}

/*! FUNCTION:  	SYSTEMIO_IsToolInstalled()
 *  SYNOPSIS:  	Checks if tool/program is installed on system.
 *                WARNING: Should only be run if tool has no side effects.
 */
bool
SYSTEMIO_IsToolInstalled( const STR tool )
{
   int ERRORCHECK_exit_code = system( tool );
   if ( ERRORCHECK_exit_code != 0 ) {
      return false;
   }
   return true;
}

/*! FUNCTION:  	SYSTEMIO_AddEnvironmentalVar()
 *  SYNOPSIS:  	Adds environmental variable 
 */
STATUS_FLAG
SYSTEMIO_AddEnvironmentalVar(    const STR   name,
                                 const STR   value )
{
   STR env_command;
   env_command    = STR_Concat( name, "=" );
   env_command    = STR_Append( env_command, value );
   putenv( env_command );
}

/*! FUNCTION:  SYSTEMIO_Wait()
 *  SYNOPSIS:  Hold for given number of <secs>.
 */
void 
SYSTEMIO_Wait( const int milli_seconds ) 
{ 
    // Converting time into milli_seconds 
    int nano_seconds =  1000 * milli_seconds; 
  
    // Storing start time 
    clock_t start_time = clock(); 
  
    // looping till required time is not achieved 
    while (clock() < start_time + nano_seconds) {}; 
} 

/*! FUNCTION:  SYSTEMIO_MakeDirectory()
 *  SYNOPSIS:  Make directory.
 *    RETURN:  Returns <STATUS_SUCCESS> on success.
 */
int 
SYSTEMIO_MakeDirectory( const char* folderpath )
{
    return STATUS_SUCCESS;
}