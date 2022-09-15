#include <stdio.h>
#include <string.h>

#include "objects/structs.h"
#include "utilities/_utilities.h"
#include "parsers/_parsers.h"
#include "pipelines/_pipelines.h"
#include "objects/_objects.h"


int main(int argc, char* argv[]) {

  if (argc < 3) {
    printf("usage: mmoreseqs <target> <query> <m8>\n");
    return EXIT_SUCCESS;
  }

  char* target_path = argv[1];
  char* query_path = argv[2];
  char* m8_path = argv[3];

  char* target_idx_path = malloc(sizeof(char) * strlen(target_path) + 4);
  char* query_idx_path = malloc(sizeof(char) * (strlen(query_path) + 4));

  strcpy(target_idx_path, target_path);
  strcpy(query_idx_path, query_path);

  strcat(target_idx_path, ".idx");
  strcat(query_idx_path, ".idx");

  WORKER* worker = WORKER_Create();
  WORKER_Init(worker);

  RNG_Init();

  COMMANDLINE_Load(worker->cmd, argc, argv);
  ARGS_SetDefaults(worker->args);

  worker->args->t_mmore_filein = target_path;
  worker->args->t_index_filein = target_idx_path;
  worker->args->q_mmore_filein = query_path;
  worker->args->q_index_filein = query_idx_path;
  worker->args->mmseqs_m8_filein = m8_path;

  worker->args->is_myout = true;
  worker->args->is_m8out = true;
  worker->args->is_myout = true;
  worker->args->is_mydom = true;
  worker->args->is_mytimeout = true;
  worker->args->is_mythreshout = true;
  worker->args->myout_fileout = "results.out";
  worker->args->m8out_fileout = "results.m8";
  worker->args->mydom_fileout = "results.domtbl";
  worker->args->mytime_fileout = "results.time";
  worker->args->mythresh_fileout = "results.thresh";

//  ARGS_Dump(worker->args, stdout);

  mmoreseqs_mmore_pipeline(worker);

  // TODO: fix ARGS_Destroy so we can still use WORKER_Destroy
//  WORKER_Destroy(worker);

  return STATUS_SUCCESS;
}

