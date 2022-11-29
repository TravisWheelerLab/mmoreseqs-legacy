#include <stdio.h>
#include <string.h>

#include "objects/structs.h"
#include "utilities/_utilities.h"
#include "parsers/_parsers.h"
#include "pipelines/_pipelines.h"
#include "objects/_objects.h"

int main(int argc, char* argv[]) {
  if (argc < 4) {
    printf("usage: mmoreseqs <query> <target> <m8> <outdir>\n");
    return EXIT_SUCCESS;
  }

  char* query_path = argv[1];
  char* target_path = argv[2];
  char* m8_path = argv[3];
  char* out_path = argv[4];

  char* myout_path = malloc(sizeof(char) * 100);
  char* m8out_path = malloc(sizeof(char) * 100);
  char* mydom_path = malloc(sizeof(char) * 100);
  char* mytime_path = malloc(sizeof(char) * 100);
  char* mythresh_path = malloc(sizeof(char) * 100);
  char* hmmerout_path = malloc(sizeof(char) * 100);

  strcat(myout_path, out_path);
  strcat(myout_path, "/results.out");
  strcat(m8out_path, out_path);
  strcat(m8out_path, "/results.m8");
  strcat(mydom_path, out_path);
  strcat(mydom_path, "/results.domtbl");
  strcat(mytime_path, out_path);
  strcat(mytime_path, "/results.time");
  strcat(mythresh_path, out_path);
  strcat(mythresh_path, "/results.thresh");
  strcat(hmmerout_path, out_path);
  strcat(hmmerout_path, "/results.hmmerout");

  printf("query: %s\n", query_path);
  printf("target: %s\n", target_path);
  printf("out: %s\n", myout_path);
  printf("m8: %s\n", m8out_path);
  printf("dom: %s\n", mydom_path);
  printf("time: %s\n", mytime_path);
  printf("thresh: %s\n", mythresh_path);
  printf("hmmerout: %s\n", hmmerout_path);
  WORKER* worker = WORKER_Create();
  WORKER_Init(worker);

  RNG_Init();

  COMMANDLINE_Load(worker->cmd, argc, argv);
  ARGS_SetDefaults(worker->args);

  worker->args->t_mmore_filein = target_path;
  worker->args->t_filein = target_path;
  worker->args->q_mmore_filein = query_path;
  worker->args->q_filein = query_path;
  worker->args->mmseqs_m8_filein = m8_path;
  worker->args->q_filetype = FILE_FASTA;
  worker->args->t_filetype = FILE_HMM;
  worker->args->is_myout = true;
  worker->args->is_m8out = true;
  worker->args->is_myout = true;
  worker->args->is_mydom = true;
  worker->args->is_mytimeout = true;
  worker->args->is_mythreshout = true;
  worker->args->is_hmmerout = true;
  worker->args->myout_fileout = myout_path;
  worker->args->m8out_fileout = m8out_path;
  worker->args->mydom_fileout = mydom_path;
  worker->args->mytime_fileout = mytime_path;
  worker->args->mythresh_fileout = mythresh_path;
  worker->args->hmmerout_fileout = hmmerout_path;


  mmoreseqs_mmore_pipeline(worker);

  // TODO: fix ARGS_Destroy so we can still use WORKER_Destroy

  //  WORKER_Destroy(worker);

  return STATUS_SUCCESS;
}
