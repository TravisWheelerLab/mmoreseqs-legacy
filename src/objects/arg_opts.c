/*******************************************************************************
 *  - FILE:  arg_opts.c
 *  - DESC:  ARG_OPTS Object.
 *           Used for Listing, Selecting, and Parsing Commandline Arguments Options.
 *  - NOTES:
 *    - WIP: Will be phased out by CLI11.
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
#include "_objects.h"

/* header */
#include "alignment.h"

/*! FUNCTION:  ARG_OPTS_Create()
 *  SYNOPSIS:  Create new ARG_OPTS object.
 */
ARG_OPTS*
ARG_OPTS_Create() {
  ARG_OPTS* arg_opts;

  arg_opts = ERROR_malloc(sizeof(ARG_OPTS));
  /* general data */
  arg_opts->opt_names = VECTOR_STR_Create();
  arg_opts->opt_desc = VECTOR_STR_Create();
  arg_opts->opt_help = VECTOR_STR_Create();
  /* option data */
  arg_opts->opt_long = VECTOR_STR_Create();
  arg_opts->opt_short = VECTOR_STR_Create();
  arg_opts->num_args = VECTOR_INT_Create();
  arg_opts->arg_index = VECTOR_INT_Create();
  /* argument data */
  arg_opts->arg_type = VECTOR_INT_Create();
  arg_opts->arg_loc = VECTOR_PTR_Create();
  arg_opts->arg_id = VECTOR_INT_Create();

  arg_opts->N_args = 0;
  arg_opts->N_opts = 0;

  VECTOR_INT_Pushback(arg_opts->arg_index, 0);

  return arg_opts;
}

/*! FUNCTION:  ARG_OPTS_Destroy()
 *  SYNOPSIS:  Destroy ARG_OPTS object.
 */
ARG_OPTS*
ARG_OPTS_Destroy(ARG_OPTS* arg_opts) {
  if (arg_opts == NULL)
    return NULL;

  /* general data */
  arg_opts->opt_names = VECTOR_STR_Destroy(arg_opts->opt_names);
  arg_opts->opt_desc = VECTOR_STR_Destroy(arg_opts->opt_desc);
  arg_opts->opt_help = VECTOR_STR_Destroy(arg_opts->opt_help);
  /* option data */
  arg_opts->opt_long = VECTOR_STR_Destroy(arg_opts->opt_long);
  arg_opts->opt_short = VECTOR_STR_Destroy(arg_opts->opt_short);
  arg_opts->num_args = VECTOR_INT_Destroy(arg_opts->num_args);
  arg_opts->arg_index = VECTOR_INT_Destroy(arg_opts->arg_index);
  /* argument data */
  arg_opts->arg_type = VECTOR_INT_Destroy(arg_opts->arg_type);
  arg_opts->arg_loc = VECTOR_PTR_Destroy(arg_opts->arg_loc);
  arg_opts->arg_id = VECTOR_INT_Destroy(arg_opts->arg_id);

  arg_opts = ERROR_free(arg_opts);

  return NULL;
}

/*! FUNCTION:  ARG_OPTS_AddOption()
 *  SYNOPSIS:  Add option to <opt_args>.
 */
void ARG_OPTS_AddOpt(ARG_OPTS* arg_opts,
                        STR name,
                        STR desc,
                        STR help,
                        STR opt_long,
                        STR opt_short,
                        int n_args,
                        PTR* arg_locs,
                        INT* arg_types) {
  int opt_id;
  int N_args;
  PTR arg_loc;
  INT arg_type;

  opt_id = VECTOR_STR_GetSize(arg_opts->opt_names);

  VECTOR_STR_Pushback(arg_opts->opt_names, name);
  VECTOR_STR_Pushback(arg_opts->opt_desc, desc);
  VECTOR_STR_Pushback(arg_opts->opt_help, help);
  VECTOR_STR_Pushback(arg_opts->opt_long, opt_long);
  VECTOR_STR_Pushback(arg_opts->opt_short, opt_short);
  VECTOR_INT_Pushback(arg_opts->num_args, n_args);

  for (int i = 0; i < n_args; i++) {
    arg_loc = arg_locs[i];
    arg_type = arg_types[i];
    VECTOR_PTR_Pushback(arg_opts->arg_loc, arg_loc);
    VECTOR_INT_Pushback(arg_opts->arg_type, arg_type);
    VECTOR_INT_Pushback(arg_opts->arg_id, opt_id);
  }

  N_args = VECTOR_STR_GetSize(arg_opts->opt_names);
  VECTOR_INT_Pushback(arg_opts->arg_index, N_args);

  arg_opts->N_opts += 1;
  arg_opts->N_args += n_args;
}

/*! FUNCTION:  ARG_OPTS_AddShOpt()
 *  SYNOPSIS:  Add shorter abbreviated option to <opt_args>.
 */
void ARG_OPTS_AddShOption(ARG_OPTS* arg_opts,
                          STR name,
                          STR opt_long,
                          int n_args,
                          PTR* arg_locs,
                          INT* arg_types) {
  int opt_id;
  int N_args;
  PTR arg_loc;
  INT arg_type;

  opt_id = VECTOR_STR_GetSize(arg_opts->opt_names);

  VECTOR_STR_Pushback(arg_opts->opt_names, name);
  VECTOR_STR_Pushback(arg_opts->opt_desc, "");
  VECTOR_STR_Pushback(arg_opts->opt_help, "");
  VECTOR_STR_Pushback(arg_opts->opt_long, opt_long);
  VECTOR_STR_Pushback(arg_opts->opt_short, "");
  VECTOR_INT_Pushback(arg_opts->num_args, n_args);

  for (int i = 0; i < n_args; i++) {
    arg_loc = arg_locs[i];
    arg_type = arg_types[i];
    VECTOR_PTR_Pushback(arg_opts->arg_loc, arg_loc);
    VECTOR_INT_Pushback(arg_opts->arg_type, arg_type);
    VECTOR_INT_Pushback(arg_opts->arg_id, opt_id);
  }

  N_args = VECTOR_STR_GetSize(arg_opts->opt_names);
  VECTOR_INT_Pushback(arg_opts->arg_index, N_args);

  arg_opts->N_opts += 1;
  arg_opts->N_args += n_args;
}


/*! FUNCTION:  ARG_OPTS_Index()
 *  SYNOPSIS:  Add option to <opt_args>.
 */
INT ARG_OPTS_Index(ARG_OPTS* arg_opts,
                   int opt_id) {
}

/*! FUNCTION:  ARG_OPTS_GetOption_byName()
 *  SYNOPSIS:  Get option <opt_name> for <arg_opts>.
 */
INT ARG_OPTS_GetOption_byName(ARG_OPTS* arg_opts,
                              STR opt_name) {
  int opt_id;
  int N_opts;

  N_opts = VECTOR_STR_GetSize(arg_opts->opt_names);
  opt_id = -1;
  for (int i = 0; i < N_opts; i++) {
    STR search_name = VECTOR_STR_Get(arg_opts->opt_names, i);
    if (strcmp(opt_name, search_name) == 0) {
      opt_id = i;
      break;
    }
  }
  return opt_id;
}

/*! FUNCTION:  ARG_OPTS_GetArgument()
 *  SYNOPSIS:  Get Argument of to <opt_args>.
 */
PTR ARG_OPTS_GetArgument(ARG_OPTS* arg_opts,
                         int arg_id) {
  PTR arg_loc;
  arg_loc = VECTOR_PTR_Get(arg_opts->arg_loc, arg_id);
  return arg_loc;
}
