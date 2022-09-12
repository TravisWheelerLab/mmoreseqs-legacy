/*******************************************************************************
 *  - FILE:  arg_opts.c
 *  - DESC:  ARG_OPTS Object.
 *******************************************************************************/

#ifndef _ARG_OPTS_H
#define _ARG_OPTS_H

/*! FUNCTION:  ARG_OPTS_Create()
 *  SYNOPSIS:  Create new ARG_OPTS object.
 */
ARG_OPTS* ARG_OPTS_Create();

/*! FUNCTION:  ARG_OPTS_Destroy()
 *  SYNOPSIS:  Destroy ARG_OPTS object.
 */
ARG_OPTS* ARG_OPTS_Destroy(ARG_OPTS* arg_opt);

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
                      INT* arg_types);

/*! FUNCTION:  ARG_OPTS_AddShOpt()
 *  SYNOPSIS:  Add shorter abbreviated option to <opt_args>.
 */
void ARG_OPTS_AddShOpt(ARG_OPTS* arg_opts,
                      STR name,
                      STR opt_long,
                      int n_args,
                      PTR* arg_locs,
                      INT* arg_types);

/*! FUNCTION:  ARG_OPTS_AddArgs()
 *  SYNOPSIS:  Get Argument of to <opt_args>.
 */
PTR ARG_OPTS_GetArgument(ARG_OPTS* arg_opts, int argid);

#endif /* ARG_OPTS_H_ */
