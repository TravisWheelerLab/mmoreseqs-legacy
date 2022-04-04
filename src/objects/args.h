/*******************************************************************************
 *  - FILE:      args.c
 *  - DESC:    ARGS Object.
 *             Used for Parsing and Storing Commandline Arguments.
 *******************************************************************************/

#ifndef _ARGS_H
#define _ARGS_H

/*! FUNCTION:  ARGS_Create()
 *  SYNOPSIS:  Create args object and allocate memory.
 */
ARGS* ARGS_Create();

/*! FUNCTION:  ARGS_Destroy()
 *  SYNOPSIS:  Destroy ARGS object and free memory.
 */
ARGS* ARGS_Destroy(ARGS* args);

#endif /* ARGS_H_ */
