/*******************************************************************************
 *  - FILE:  arg_parser.h
 *  - DESC:  Parses command line arguments.
 *******************************************************************************/

#ifndef _ARG_PARSER
#define _ARG_PARSER

/*! FUNCTION:  ARGS_Parse()
 *  SYNOPSIS:  Parses Arguments from the command line.
 */
void ARGS_Parse(ARGS* args,
                int argc,
                char* argv[],
                COMMANDLINE* cmd,
                ARG_OPTS* arg_opts);

/*! FUNCTION:  ARGS_SetDefaults()
 *  SYNOPSIS:  Set default arguments.
 */
void ARGS_SetDefaults(ARGS* args);

/*! FUNCTION:  ARGS_SetDefaults()
 *  SYNOPSIS:  Output arguments to <fp>.
 */
void ARGS_Dump(ARGS* args, FILE* fp);

/*! FUNCTION:  ARGS_FindFiletype()
 *  SYNOPSIS:  Examines target and query, and finds the type of the files (by
 * extension).
 */
FILETYPE
ARGS_FindFiletype(STR filename);

/*! FUNCTION:  ARGS_HelpInfo()
 *  SYNOPSIS:  Output help info.
 */
void ARGS_HelpInfo();

/*! FUNCTION:  ARGS_CommandHelpInfo()
 *  SYNOPSIS:  Output command help info.
 */
void ARGS_CommandHelpInfo(ARGS* args);

/*! FUNCTION:  ARGS_VersionInfo()
 *  SYNOPSIS:  Output version info.
 */
void ARGS_VersionInfo();

/*! FUNCTION:  ARGS_Default_Opts()
 *  SYNOPSIS:  Output Default Options.
 */
STATUS_FLAG
ARGS_ParseCommand(ARGS* args,
                  int argc_p, 
                  char* argv[],
                  int* arg_cur_p);

/*! FUNCTION:  ARGS_Default_Opts()
 *  SYNOPSIS:  Output Default Options.
 */
STATUS_FLAG
ARGS_ParseOptions(ARGS* args,
                  int argc, 
                  char* argv[],
                  int* arg_cur_p);

/*! FUNCTION:  ARGS_Default_Opts()
 *  SYNOPSIS:  Output Default Options.
 */
STATUS_FLAG ARGS_ParseOptionsArgp(ARGS* args,
                                  int argc,
                                  char* argv[]);

#endif /* _ARG_PARSER */
