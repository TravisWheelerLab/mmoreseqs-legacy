/*******************************************************************************
 *  - FILE:      args.c
 *  - DESC:    ARGS Object.
 *             Used for Parsing and Storing Commandline Arguments.
 *******************************************************************************/

#ifndef _ARGS_H
#define _ARGS_H

/*! FUNCTION:  ARGS_Destroy()
 *  SYNOPSIS:  Destroy ARGS object and free memory.
 */
ARGS* ARGS_Destroy(ARGS* args);

#endif  //_ARGS_H