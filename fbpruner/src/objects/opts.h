/*******************************************************************************
 *  FILE:      opts.c
 *  PURPOSE:   OPTS Object. Used for Parsing Commandline Args.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _OPTS_H
#define _OPTS_H

/*
 *  FUNCTION:  OPTS_Create()
 *  SYNOPSIS:  
 */
OPTS* 
OPTS_Create();

/*
 *  FUNCTION:  OPTS_Destroy()
 *  SYNOPSIS:  
 */
OPTS* 
OPTS_Destroy( OPTS* opts );

#endif /* _OPTS_H */