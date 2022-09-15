/*******************************************************************************
 *  - FILE:      commandline.c
 *  - DESC:    COMMANDLINE Object.
 *             Stores commandline arguments.
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
COMMANDLINE_Create() {
  COMMANDLINE* cmd;
  cmd = ERROR_malloc(sizeof(COMMANDLINE));

  cmd->cur_cmd = 0;
  cmd->cur_opt = 0;
  cmd->cur_arg = 0;

  cmd->cmds = VECTOR_STR_Create();
  cmd->options = VECTOR_STR_Create();
  cmd->arguments = VECTOR_STR_Create();
  cmd->argument_opts = VECTOR_INT_Create();
  cmd->argument_index = VECTOR_INT_Create();

  return cmd;
}

/*! FUNCTION:  COMMANDLINE_Destroy()
 *  SYNOPSIS:  Destroy commandline object.
 */
COMMANDLINE*
COMMANDLINE_Destroy(COMMANDLINE* cmd) {
  if (cmd == NULL)
    return NULL;

  cmd->cmds = VECTOR_STR_Destroy(cmd->cmds);
  cmd->options = VECTOR_STR_Destroy(cmd->options);
  cmd->arguments = VECTOR_STR_Destroy(cmd->arguments);
  cmd->argument_opts = VECTOR_INT_Destroy(cmd->argument_opts);
  cmd->argument_index = VECTOR_INT_Destroy(cmd->argument_index);
  cmd = ERROR_free(cmd);

  return NULL;
}

/*! FUNCTION:  COMMANDLINE_Reuse()
 *  SYNOPSIS:  Reuse commandline object.
 */
void COMMANDLINE_Reuse(COMMANDLINE* cmd) {
  cmd->cur_cmd = 0;
  cmd->cur_opt = 0;
  cmd->cur_arg = 0;

  VECTOR_STR_Reuse(cmd->cmds);
  VECTOR_STR_Reuse(cmd->options);
  VECTOR_STR_Reuse(cmd->arguments);
  VECTOR_INT_Reuse(cmd->argument_opts);
  VECTOR_INT_Reuse(cmd->argument_index);
}

/*! FUNCTION:  COMMANDLINE_Get()
 *  SYNOPSIS:  Get <i>th commandline entry from <cmd>.
 *             Caller must load
 */
STR COMMANDLINE_Get(COMMANDLINE* cmd,
                    int i) {
  return VECTOR_STR_Get(cmd->cmds, i);
}

/*! FUNCTION:  COMMANDLINE_Load()
 *  SYNOPSIS:  Load commandline into <cmd> object.
 */
void COMMANDLINE_Load(COMMANDLINE* cmd,
                      int argc,
                      STR argv[]) {
  int long_flag;
  int short_flag;
  int opt_num = 0;
  int arg_num = 0;

  VECTOR_INT_Pushback(cmd->argument_index, 0);
  VECTOR_STR_Pushback(cmd->options, NULL);

  /* all commandline args */
  for (int i = 0; i < argc; i++) {
    char* arg = argv[i];

    /* add to full list */
    VECTOR_STR_Pushback(cmd->cmds, arg);

    /* check if long flag option */
    if (long_flag = STR_StartsWith(arg, "--"), long_flag == 0) {
      // printf("long_flag: %d, %s\n", long_flag, arg );
      opt_num += 1;
      VECTOR_STR_Pushback(cmd->options, arg);
      VECTOR_INT_Pushback(cmd->argument_index, arg_num);
    }
    /* check if short flag option */
    elif (short_flag = STR_StartsWith(arg, "-"), short_flag == 0) {
      // printf("short_flag: %d, %s\n", short_flag, arg );
      opt_num += 1;
      VECTOR_STR_Pushback(cmd->options, arg);
      VECTOR_INT_Pushback(cmd->argument_index, arg_num);
    }
    /* else it is an argument */
    else {
      // printf("argument: %s\n", arg );
      arg_num += 1;
      VECTOR_STR_Pushback(cmd->arguments, arg);
      VECTOR_INT_Pushback(cmd->argument_opts, opt_num);
    }
  }
  VECTOR_INT_Pushback(cmd->argument_index, arg_num);
}

/*! FUNCTION:  COMMANDLINE_GetNumOpts()
 *  SYNOPSIS:  Get number of options.
 */
size_t
COMMANDLINE_GetNumCmds(COMMANDLINE* cmd) {
  size_t N;
  N = VECTOR_STR_GetSize(cmd->cmds);
  return N;
}


/*! FUNCTION:  COMMANDLINE_SimpleDump()
 *  SYNOPSIS:  Output <cmd> to <fp>.
 */
void COMMANDLINE_SimpleDump(COMMANDLINE* cmd,
                            FILE* fp) {
  int N_cmds;
  STR my_cmd;

  N_cmds = COMMANDLINE_GetNumCmds(cmd);
  fprintf(fp, "# === COMMAND: ");

  for (int i = 0; i < N_cmds; i++) {
    my_cmd = COMMANDLINE_Get(cmd, i);
    printf("%s ", my_cmd);
  }
  printf("\n");
}
