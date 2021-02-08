/*******************************************************************************
 *  FILE:      arg_parser.h
 *  PURPOSE:   Parses command line arguments.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _ARG_PARSER
#define _ARG_PARSER

/* parses arguments from the command line */
void ARGS_Parse(  ARGS*   args,
				      int     argc, 
                  char*   argv[] );

/* set default arguments */
void ARGS_SetDefaults( ARGS* args );

/* sends ARGS data to FILE POINTER */
void ARGS_Dump( ARGS*    args,
                FILE*    fp );

/* examines target and query, and finds the type of the files */
int ARGS_Find_FileType( char* _filename_ );

/* output help info */
void ARGS_Help_Info();

#endif /* _ARG_PARSER */