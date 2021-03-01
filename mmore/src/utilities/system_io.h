/*******************************************************************************
 *  FILE:      system_io.h
 *  PURPOSE:   Tools for file input/output.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:
 *******************************************************************************/

#ifndef _SYSTEMIO_H
#define _SYSTEMIO_H

#include "../objects/structs.h"

/*! FUNCTION:  	SYSTEMIO_Exists()
 *  SYNOPSIS:  	Checks whether <filename> exists.
 */
bool
SYSTEMIO_Exists( const STR filename );

/*! FUNCTION:  	SYSTEMIO_HasWritePermission()
 *  SYNOPSIS:  	Checks whether user has permission to write to <filename>.
 */
bool
SYSTEMIO_HasWritePermission( const STR filename );

/*! FUNCTION:  	SYSTEMIO_HasReadPermission()
 *  SYNOPSIS:  	Returns whether user has permission to read from <filename>.
 */
bool
SYSTEMIO_HasReadPermission( const STR filename );

/*! FUNCTION:  	SYSTEMIO_GetDirectory()
 *  SYNOPSIS:  	Returns current working directory.
 *                Malloc's string data if passed NULL pointer.
 */
STR 
SYSTEMIO_GetDirectory( STR old_str );

/*! FUNCTION:  	SYSTEMIO_IsToolInstalled()
 *  SYNOPSIS:  	Checks if tool/program is installed on system.
 *                WARNING: Should only be run if tool has no side effects.
 */
bool
SYSTEMIO_IsToolInstalled( const STR tool );

/*! FUNCTION:  	SYSTEMIO_AddEnvironmentalVar()
 *  SYNOPSIS:  	Adds environmental variable.
 */
STATUS_FLAG
SYSTEMIO_AddEnvironmentalVar(    const STR   name,
                                 const STR   value );

/*! FUNCTION:  SYSTEMIO_Wait()
 *  SYNOPSIS:  Hold for given number of <secs>.
 */
void 
SYSTEMIO_Wait( int milli_seconds );

/*! FUNCTION:  SYSTEMIO_MakeDirectory()
 *  SYNOPSIS:  Make directory.
 *    RETURN:  Returns <STATUS_SUCCESS> on success.
 */
int 
SYSTEMIO_MakeDirectory( const char* folderpath );

#endif /* _SYSTEMIO_H */