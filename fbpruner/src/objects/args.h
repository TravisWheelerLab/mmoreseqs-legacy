/*******************************************************************************
 *  FILE:      args.c
 *  PURPOSE:   ARGS Object. Used for Parsing Commandline Args.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _ARGS_H
#define _ARGS_H

/*
 *  FUNCTION:  ARGS_Create()
 *  SYNOPSIS:  
 */
ARGS* 
ARGS_Create();

/*
 *  FUNCTION:  ARGS_Destroy()
 *  SYNOPSIS:  
 */
void* 
ARGS_Destroy(ARGS* args);

#endif /* ARGS_H_ */