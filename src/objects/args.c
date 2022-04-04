/*******************************************************************************
 *  - FILE:      args.c
 *  - DESC:    ARGS Object. Used for Parsing Commandline Args.
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

/*! FUNCTION:  ARGS_Create()
 *  SYNOPSIS:
 */
ARGS* ARGS_Create() {
  ARGS* args = NULL;
  args = ERROR_malloc(sizeof(ARGS));

  /* commandline */
  args->cmdline = NULL;
  args->opts = NULL;
  /* pipeline */
  args->pipeline_mode = -1;
  args->pipeline_name = NULL;
  args->search_name = NULL;
  /* folders */
  args->tmp_folderpath = NULL;
  args->dbg_folderpath = NULL;
  args->prep_folderpath = NULL;
  /* input */
  args->t_filein = NULL;
  args->q_filein = NULL;
  /* mmore input */
  args->t_mmore_filein = NULL;
  args->t_mmore_p_filein = NULL;
  args->t_mmore_s_filein = NULL;
  args->q_mmore_filein = NULL;
  /* mmseqs input */
  args->t_mmseqs_filein = NULL;
  args->t_mmseqs_p_filein = NULL;
  args->t_mmseqs_s_filein = NULL;
  args->q_mmseqs_filein = NULL;
  /* mmseqs output */
  args->mmseqs_m8out = NULL;
  /* index input */
  args->t_index_filein = NULL;
  args->q_index_filein = NULL;
  /* results input */
  args->mmseqs_m8_filein = NULL;
  args->hitlist_filein = NULL;
  /* output */
  args->stdout_fileout = NULL;
  args->stderr_fileout = NULL;
  args->allout_fileout = NULL;
  args->hmmerout_fileout = NULL;
  args->m8out_fileout = NULL;
  args->myout_fileout = NULL;
  args->mydom_fileout = NULL;
  args->mytime_fileout = NULL;
  args->mythresh_fileout = NULL;
  args->customout_fileout = NULL;

  return args;
}

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
