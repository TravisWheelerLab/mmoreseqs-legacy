/*******************************************************************************
 *  - FILE:      commandline.h
 *  - DESC:    COMMANDLINE Object.
 *             Stores commandline arguments.
 *******************************************************************************/

#ifndef _COMMANDLINE_H
#define _COMMANDLINE_H

/*! FUNCTION:  COMMANDLINE_Create()
 *  SYNOPSIS:  Create new commandline object.
 */
COMMANDLINE* COMMANDLINE_Create();

/*! FUNCTION:  COMMANDLINE_Destroy()
 *  SYNOPSIS:  Destroy commandline object.
 */
COMMANDLINE* COMMANDLINE_Destroy(COMMANDLINE* cmd);

/*! FUNCTION:  COMMANDLINE_Reuse()
 *  SYNOPSIS:  Reuse commandline object.
 */
void COMMANDLINE_Reuse(COMMANDLINE* cmd);

/*! FUNCTION:  COMMANDLINE_LoadArgs()
 *  SYNOPSIS:  Load commandline into <cmd> object.
 */
void COMMANDLINE_Load(COMMANDLINE* cmd, int argc, STR argv[]);

/*! FUNCTION:  COMMANDLINE_GetNumOpts()
 *  SYNOPSIS:  Get number of total options.
 */
size_t COMMANDLINE_GetNumOpts(COMMANDLINE* cmd);

/*! FUNCTION:  COMMANDLINE_GetNumArgs()
 *  SYNOPSIS:  Get number of total arguments.
 */
size_t COMMANDLINE_GetNumArgs(COMMANDLINE* cmd);

/*! FUNCTION:  COMMANDLINE_GetNumArgs()
 *  SYNOPSIS:  Get argument range for <opt_id> option.
 */
RANGE
COMMANDLINE_GetRangeArgs_byOpt(COMMANDLINE* cmd, int opt_id);

/*! FUNCTION:  COMMANDLINE_Dump()
 *  SYNOPSIS:  Output <cmd> to <fp>.
 */
void COMMANDLINE_Dump(COMMANDLINE* cmd, FILE* fp);

/*! FUNCTION:  COMMANDLINE_SimpleDump()
 *  SYNOPSIS:  Output <cmd> to <fp>.
 */
void COMMANDLINE_SimpleDump(COMMANDLINE* cmd, FILE* fp);

#endif /* COMMANDLINE_H_ */
