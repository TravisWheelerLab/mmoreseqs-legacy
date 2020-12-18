/*******************************************************************************
 *  FILE:      arg_parser.h
 *  PURPOSE:   Parses command line arguments.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:
 *******************************************************************************/


/* parses arguments from the command line */
void ARGS_Parse(  ARGS*   args,
				      int     argc, 
                  char*   argv[] );

/* set default arguments */
void ARGS_Set_Defaults( ARGS* args );

/* sends ARGS data to FILE POINTER */
void ARGS_Dump( ARGS*    args,
                FILE*    fp );

/* examines target and query, and finds the type of the files */
int ARGS_Find_FileType( char* _filename_ );

/* output help info */
void ARGS_Help_Info();