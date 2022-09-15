/*******************************************************************************
 *  - FILE:      args.c
 *  - DESC:    ARGS Object. Used for Parsing Commandline Args.
 *******************************************************************************/

/* imports */
#include <string.h>

/* local imports */
#include "structs.h"
#include "../utilities/_utilities.h"
#include "_objects.h"

/*! FUNCTION:  ARGS_Destroy()
 *  SYNOPSIS:
 */
ARGS* ARGS_Destroy(ARGS* args) {
  if (args == NULL)
    return args;

  /* commandline */
  STR_Destroy(args->cmdline);
  STR_Destroy(args->opts);
  /* program */
  STR_Destroy(args->mmoreseqs_program);
  STR_Destroy(args->mmseqs_program);
  STR_Destroy(args->hmmer_program);
  STR_Destroy(args->mmoreseqs_scripts);
  /* pipeline */
  STR_Destroy(args->pipeline_name);
  STR_Destroy(args->search_name);
  /* folders */
  STR_Destroy(args->tmp_folderpath);
  STR_Destroy(args->dbg_folderpath);
  STR_Destroy(args->prep_folderpath);
  /* input file */
  STR_Destroy(args->t_filein);
  STR_Destroy(args->q_filein);
  /* mmore input */
  STR_Destroy(args->t_mmore_filein);
  STR_Destroy(args->t_mmore_p_filein);
  STR_Destroy(args->t_mmore_s_filein);
  STR_Destroy(args->q_mmore_filein);
  /* mmseqs input */
  STR_Destroy(args->t_mmseqs_filein);
  STR_Destroy(args->t_mmseqs_p_filein);
  STR_Destroy(args->t_mmseqs_s_filein);
  STR_Destroy(args->q_mmseqs_filein);
  /* mmseqs output */
  STR_Destroy(args->mmseqs_m8out);
  /* index input */
  STR_Destroy(args->t_index_filein);
  STR_Destroy(args->q_index_filein);
  /* results input */
  STR_Destroy(args->mmseqs_m8_filein);
  STR_Destroy(args->hitlist_filein);
  /* prep-able files */
  STR_Destroy(args->target_prep);
  STR_Destroy(args->query_prep);
  STR_Destroy(args->link_target_mmore_prep);
  STR_Destroy(args->link_query_mmore_prep);
  STR_Destroy(args->link_target_mmseqs_prep);
  STR_Destroy(args->link_query_mmseqs_prep);
  /* output file */
  STR_Destroy(args->stdout_fileout);
  STR_Destroy(args->stderr_fileout);
  STR_Destroy(args->allout_fileout);
  STR_Destroy(args->hmmerout_fileout);
  STR_Destroy(args->m8out_fileout);
  STR_Destroy(args->myout_fileout);
  STR_Destroy(args->mydom_fileout);
  STR_Destroy(args->mytime_fileout);
  STR_Destroy(args->mythresh_fileout);
  STR_Destroy(args->customout_fileout);

  args = ERROR_free(args);
  return args;
}
