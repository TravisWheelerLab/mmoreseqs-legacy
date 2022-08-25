/*******************************************************************************
 *  - FILE:  application.c
 *  - DESC:  Entry Point to Application: Worker/Debugger Initialization, Argument Parsing, Chooses Pipeline
 *******************************************************************************/

#ifndef _APPLICATION_H
#define _APPLICATION_H

#include "objects/structs.h"

/*! FUNCTION:  	run_application()
 *  SYNOPSIS:  	Main entry point into application.
 */
STATUS_FLAG
run_application(int argc, char* argv[]);

#endif /* _APPLICATION_H */
