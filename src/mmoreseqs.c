/*******************************************************************************
 *  - FILE:  mmoreseqs.c
 *  - DESC:  Entry Point to Application: Worker/Debugger Initialization, Argument Parsing, Chooses Pipeline
 *******************************************************************************/

/* import stdlib */
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
// #include <dos.h> /* date function */

/* import local libraries */
#include "easel.h"

/* include local files */
#include "objects/structs.h"
#include "utilities/_utilities.h"
#include "parsers/_parsers.h"
#include "pipelines/_pipelines.h"
#include "objects/_objects.h"

/* header */
#include "application.h"

/* === MAIN ENTRY-POINT TO PROGRAM === */
STATUS_FLAG
main(int argc, char* argv[]) {
  STATUS_FLAG return_status = APPLICATION_Run(argc, argv);
  ERRORCHECK_exit(return_status);
}
