/*******************************************************************************
 *  FILE:      main.c
 *  PURPOSE:   Main Method, Parses Command Line Arguments, then
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/



/* === FUNCTIONS === */

/* Parses Arguments from the command line */
ARGS* ARGS_Parse(int     argc, 
                 char*   argv[]);

/* SET DEFAULT ARGUMENTS (for testing) */
void  ARGS_Set_Defaults(ARGS* args);

/* sends ARGS data to FILE POINTER */
void ARGS_Dump(ARGS*    args,
               FILE*    fp);

/* examines target and query, and finds the type of the files */
int determine_FileType( char* _filename_ );