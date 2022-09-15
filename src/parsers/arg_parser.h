/*******************************************************************************
 *  - FILE:  arg_parser.h
 *  - DESC:  Parses command line arguments.
 *******************************************************************************/

#ifndef _ARG_PARSER
#define _ARG_PARSER


/*! FUNCTION:  ARGS_SetDefaults()
 *  SYNOPSIS:  Set default arguments.
 */
void ARGS_SetDefaults(ARGS* args);

/*! FUNCTION:  ARGS_SetDefaults()
 *  SYNOPSIS:  Output arguments to <fp>.
 */
void ARGS_Dump(ARGS* args, FILE* fp);



#endif /* _ARG_PARSER */
