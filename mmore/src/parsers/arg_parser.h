/*******************************************************************************
 *  FILE:      arg_parser.h
 *  PURPOSE:   Parses command line arguments.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _ARG_PARSER
#define _ARG_PARSER

/*! FUNCTION:  ARGS_Parse()
 *  SYNOPSIS:  Parses Arguments from the command line.
 */
void   
ARGS_Parse(    ARGS*          args,
               int            argc, 
               char*          argv[],
               COMMANDLINE*   cmd,
               ARG_OPTS*      arg_opts );

/*! FUNCTION:  ARGS_SetDefaults()
 *  SYNOPSIS:  Set default arguments.
 */
void  
ARGS_SetDefaults( ARGS*    args );

/*! FUNCTION:  ARGS_MainArg_Parser()
 *  SYNOPSIS:  Set default arguments.
 */
void
ARGS_MainArg_Parser(    ARGS*          args,
                        COMMANDLINE*   cmd );

/*! FUNCTION:  ARGS_ArgOpt_Parser()
 *  SYNOPSIS:  Set default arguments.
 */
ARGS_ArgOpt_Parser(  ARGS*          args, 
                     COMMANDLINE*   cmd, 
                     STR            opt,
                     DATATYPE       dtype );

/*! FUNCTION:  ARGS_SetDefaults()
 *  SYNOPSIS:  Output arguments to <fp>.
 */
void 
ARGS_Dump(     ARGS*    args,
               FILE*    fp );

/*! FUNCTION:  ARGS_Find_Filetype()
 *  SYNOPSIS:  Examines target and query, and finds the type of the files (by extension).
 */
FILETYPE 
ARGS_FindFiletype( STR    filename );

/*! FUNCTION:  ARGS_Help_Info()
 *  SYNOPSIS:  Output help info.
 */
void 
ARGS_Help_Info();

/*! FUNCTION:  ARGS_Version_Info()
 *  SYNOPSIS:  Output version info.
 */
void 
ARGS_Version_Info();

/*! FUNCTION:  ARGS_Default_Opts()
 *  SYNOPSIS:  Output Default Options.
 */
STATUS_FLAG
ARGS_SetOptions(  ARGS*       args,
                  ARG_OPTS*   arg_opts );

#endif /* _ARG_PARSER */