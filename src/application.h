/*******************************************************************************
 *  - FILE:  application.c
 *  - DESC:  Entry Point to Application: Worker/Debugger Initialization, Argument Parsing, Chooses Pipeline
 *******************************************************************************/

#ifndef _APPLICATION_H
#define _APPLICATION_H

#include "objects/structs.h"

/*! FUNCTION:  	APPLICATION_Run()
 *  SYNOPSIS:  	Main entry point into application.
 */
STATUS_FLAG
APPLICATION_Run(int argc, char* argv[]);

#endif /* _APPLICATION_H */
