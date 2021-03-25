/*******************************************************************************
 *  FILE:      commandline.c
 *  PURPOSE:   COMMANDLINE Object. 
 *             Stores commandline arguments.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/* local imports */
#include "structs.h"
#include "../utilities/_utilities.h"

/* header */
#include "_objects.h"
#include "commandline.h"

/*! FUNCTION:  COMMANDLINE_Create()
 *  SYNOPSIS:  Create new commandline object.
 */
COMMANDLINE* 
COMMANDLINE_Create( )
{
   COMMANDLINE* cmd;
   cmd = ERROR_malloc( sizeof(COMMANDLINE) );

   cmd->cur_cmd         = 0;
   cmd->cur_opt         = 0;
   cmd->cur_arg         = 0;

   cmd->cmds            = VECTOR_STR_Create();
   cmd->options         = VECTOR_STR_Create();
   cmd->arguments       = VECTOR_STR_Create();
   cmd->argument_opts   = VECTOR_INT_Create();
   cmd->argument_index  = VECTOR_INT_Create();

   return cmd;
}

/*! FUNCTION:  COMMANDLINE_Destroy()
 *  SYNOPSIS:  Destroy commandline object.
 */
COMMANDLINE* 
COMMANDLINE_Destroy( COMMANDLINE* cmd )
{
   if ( cmd == NULL ) return NULL;

   cmd->cmds            = VECTOR_STR_Destroy( cmd->cmds );
   cmd->options         = VECTOR_STR_Destroy( cmd->options );
   cmd->arguments       = VECTOR_STR_Destroy( cmd->arguments );
   cmd->argument_opts   = VECTOR_INT_Destroy( cmd->argument_opts );
   cmd->argument_index  = VECTOR_INT_Destroy( cmd->argument_index );
   cmd                  = ERROR_free( cmd );

   return NULL;
}

/*! FUNCTION:  COMMANDLINE_Reuse()
 *  SYNOPSIS:  Reuse commandline object.
 */
void 
COMMANDLINE_Reuse( COMMANDLINE* cmd )
{
   cmd->cur_cmd = 0;
   cmd->cur_opt = 0;
   cmd->cur_arg = 0;

   VECTOR_STR_Reuse( cmd->cmds );
   VECTOR_STR_Reuse( cmd->options );
   VECTOR_STR_Reuse( cmd->arguments );
   VECTOR_INT_Reuse( cmd->argument_opts );
   VECTOR_INT_Reuse( cmd->argument_index );
}

/*! FUNCTION:  COMMANDLINE_Get()
 *  SYNOPSIS:  Get <i>th commandline entry from <cmd>.
 *             Caller must load
 */
STR 
COMMANDLINE_Get(  COMMANDLINE*   cmd,
                  int            i )
{
   return VECTOR_STR_Get( cmd->cmds, i );
}

/*! FUNCTION:  COMMANDLINE_Load()
 *  SYNOPSIS:  Load commandline into <cmd> object.
 */
void 
COMMANDLINE_Load(    COMMANDLINE*   cmd,
                     int            argc,
                     STR            argv[] )
{
   int   long_flag;
   int   short_flag;
   int   opt_num        = 0;
   int   arg_num        = 0;

   VECTOR_INT_Pushback( cmd->argument_index, 0 );
   VECTOR_STR_Pushback( cmd->options, NULL );

   /* all commandline args */
   for (int i = 0; i < argc; i++) 
   {
      char* arg = argv[i];

      /* add to full list */
      VECTOR_STR_Pushback( cmd->cmds, arg );

      /* check if long flag option */
      if ( long_flag = STR_StartsWith( arg, "--" ), long_flag == 0 ) {
         // printf("long_flag: %d, %s\n", long_flag, arg );
         opt_num += 1;
         VECTOR_STR_Pushback( cmd->options, arg );
         VECTOR_INT_Pushback( cmd->argument_index, arg_num );
      }
      /* check if short flag option */
      elif ( short_flag = STR_StartsWith( arg, "-" ), short_flag == 0 ) {
         // printf("short_flag: %d, %s\n", short_flag, arg );
         opt_num += 1;
         VECTOR_STR_Pushback( cmd->options, arg );
         VECTOR_INT_Pushback( cmd->argument_index, arg_num );
      }
      /* else it is an argument */
      else {
         // printf("argument: %s\n", arg );
         arg_num += 1;
         VECTOR_STR_Pushback( cmd->arguments, arg );
         VECTOR_INT_Pushback( cmd->argument_opts, opt_num );
      }
   }
   VECTOR_INT_Pushback( cmd->argument_index, arg_num );
}

/*! FUNCTION:  COMMANDLINE_GetNumOpts()
 *  SYNOPSIS:  Get number of options.
 */
size_t 
COMMANDLINE_GetNumCmds( COMMANDLINE*   cmd )
{
   size_t N;
   N = VECTOR_STR_GetSize( cmd->cmds );
   return N;
}

/*! FUNCTION:  COMMANDLINE_GetNumOpts()
 *  SYNOPSIS:  Get number of options.
 */
size_t 
COMMANDLINE_GetNumOpts( COMMANDLINE*   cmd )
{
   size_t N;
   N = VECTOR_STR_GetSize( cmd->options );
   return N;
}

/*! /*! FUNCTION:  COMMANDLINE_GetNumArgs()
 *  SYNOPSIS:  Get number of arguments.
 */
size_t 
COMMANDLINE_GetNumArgs( COMMANDLINE*   cmd )
{
   size_t N;
   N = VECTOR_STR_GetSize( cmd->arguments );
   return N;
}

/*! FUNCTION:  COMMANDLINE_GetRangeArgs_byOpt()
 *  SYNOPSIS:  Get argument range for <opt_id> option.
 */
RANGE 
COMMANDLINE_GetRangeArgs_byOpt(  COMMANDLINE*   cmd,
                                 int            opt_id )
{
   RANGE rng;
   rng.beg  = VECTOR_INT_Get( cmd->argument_index, opt_id );
   rng.end  = VECTOR_INT_Get( cmd->argument_index, opt_id+1 );
   return rng;
}

/*! FUNCTION:  COMMANDLINE_GetNumArgs_byOpt()
 *  SYNOPSIS:  Get argument range for <opt_id> option.
 */
int 
COMMANDLINE_GetNumArgs_byOpt(    COMMANDLINE*   cmd,
                                 int            opt_id )
{
   RANGE    rng_args;
   int      num_args;
   rng_args    = COMMANDLINE_GetRangeArgs_byOpt( cmd, opt_id );
   num_args    = rng_args.end - rng_args.beg;
   return num_args;  
}

/*! FUNCTION:  COMMANDLINE_Dump()
 *  SYNOPSIS:  Output <cmd> to <fp>.
 */
void 
COMMANDLINE_Dump( COMMANDLINE*  cmd,
                  FILE*         fp )
{
   int      N_opts, N_args;
   STR      opt;
   STR      arg;
   RANGE    arg_rng; 
   int      i,j;
   int      j_arg;
   
   N_opts = COMMANDLINE_GetNumOpts( cmd );
   N_args = COMMANDLINE_GetNumArgs( cmd );

   fprintf(fp, "=== COMMANDLINE ===:\n");
   VECTOR_STR_Dump_byOpt( cmd->options,         ",",  "OPTIONS",        stdout );
   VECTOR_STR_Dump_byOpt( cmd->arguments,       ",",  "ARGUMENTS",      stdout );
   VECTOR_INT_Dump_byOpt( cmd->argument_index,  ",",  "ARG_INDEX",      stdout );
   VECTOR_INT_Dump_byOpt( cmd->argument_opts,   ",",  "ARG_OPTION_ID",  stdout );

   for (i = 0; i < N_opts; i++) 
   {
      opt         = VECTOR_STR_Get( cmd->options, i );
      arg_rng.beg = VECTOR_INT_Get( cmd->argument_index, i );  
      arg_rng.end = VECTOR_INT_Get( cmd->argument_index, i+1 );

      fprintf(fp, "OPT[%d]: %s", i, opt );
      fprintf(fp, "\t{ ");

      for (j = arg_rng.beg; j < arg_rng.end; j++)
      {
         arg      = VECTOR_STR_Get( cmd->arguments, j );
         j_arg    = j - arg_rng.beg;

         fprintf(fp, "[%d]: %s", j_arg, arg );
         if ( j < arg_rng.end - 1 ) {
            fprintf(fp, ", ");
         }
      }  
      fprintf(fp, " }\n");
   }
}

/*! FUNCTION:  COMMANDLINE_SimpleDump()
 *  SYNOPSIS:  Output <cmd> to <fp>.
 */
void 
COMMANDLINE_SimpleDump( COMMANDLINE*  cmd,
                        FILE*         fp )
{
   int      N_cmds;
   STR      my_cmd;

   N_cmds = COMMANDLINE_GetNumCmds( cmd );
   fprintf(fp, "# COMMAND: ");

   for (int i = 0; i < N_cmds; i++) 
   {
      my_cmd = COMMANDLINE_Get( cmd, i );
      printf("%s ", my_cmd);
   }
   printf("\n");
}