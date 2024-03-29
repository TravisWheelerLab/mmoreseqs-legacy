/*******************************************************************************
 *  - FILE:      hmm_profile.c
 *  - DESC:    HMM_PROFILE Object.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

/* local imports */
#include "structs.h"
#include "../utilities/_utilities.h"
#include "_objects.h"

/* header */
#include "hmm_profile.h"

/**  FUNCTION:  HMM_PROFILE_Create()
 *   SYNOPSIS:  Constructor for HMM_PROFILE.
 *   RETURN:    Pointer to HMM_PROFILE object.
 */
HMM_PROFILE*
HMM_PROFILE_Create() {
  HMM_PROFILE* prof = NULL;

  prof = (HMM_PROFILE*)ERROR_malloc(sizeof(HMM_PROFILE));

  prof->filepath = NULL;
  prof->format = NULL;
  prof->name = NULL;
  prof->acc = NULL;
  prof->desc = NULL;
  prof->alph = NULL;

  prof->consensus = VECTOR_CHAR_Create();
  prof->is_consensus = false;

  prof->msv_dist = (DIST_PARAM){0.0, 0.0};
  prof->viterbi_dist = (DIST_PARAM){0.0, 0.0};
  prof->forward_dist = (DIST_PARAM){0.0, 0.0};

  prof->N = 0;
  prof->Nalloc = 0;
  prof->alph_leng = 20;

  prof->bg_model = NULL;
  prof->hmm_model = NULL;

  /* temporary store of full sequence */
  prof->N_full = -1;
  prof->hmm_model_full = NULL;

  prof->bg_model = HMM_COMPO_Create();

  return prof;
}

/**  FUNCTION:  HMM_PROFILE_Destroy()
 *   SYNOPSIS:  Destructor for HMM_PROFILE. Frees all associated memory.
 *   RETURN:    NULL pointer.
 */
HMM_PROFILE*
HMM_PROFILE_Destroy(HMM_PROFILE* prof) {
  if (prof == NULL)
    return prof;

  prof->filepath = ERROR_free(prof->filepath);
  prof->format = ERROR_free(prof->format);
  prof->name = ERROR_free(prof->name);
  prof->acc = ERROR_free(prof->acc);
  prof->desc = ERROR_free(prof->desc);
  prof->alph = ERROR_free(prof->alph);

  VECTOR_CHAR_Destroy(prof->consensus);

  ERROR_free(prof->bg_model);
  ERROR_free(prof->hmm_model);

  prof = ERROR_free(prof);
  return prof;
}

/* reuse profile by setting length of length of profile to zero */
/**  FUNCTION:  HMM_PROFILE_Reuse()
 *   SYNOPSIS:
 */
int HMM_PROFILE_Reuse(HMM_PROFILE* prof) {
  prof->N = 0;
  VECTOR_CHAR_Reuse(prof->consensus);
  prof->is_consensus = false;

  return STATUS_SUCCESS;
}

/* create backround hmm composition from hardcoded background frequencies */
/**  FUNCTION:  HMM_PROFILE_Reuse()
 *   SYNOPSIS:
 */
HMM_COMPO*
HMM_COMPO_Create() {
  HMM_COMPO* bg = ERROR_malloc(sizeof(HMM_COMPO));

  /* set background frequencies from BG_MODEL (imported from HMMER) */
  for (int i = 0; i < NUM_AMINO; i++) {
    bg->freq[i] = BG_MODEL[i];
  }

  return bg;
}

/* Set text to given HMM_PROFILE field */
/**  FUNCTION:  HMM_PROFILE_Reuse()
 *   SYNOPSIS:
 */
void HMM_PROFILE_SetTextField(char** prof_field,
                              char* text) {
  *prof_field = STR_Set(*prof_field, text);
}

/* Set HMM Model Length and allocate memory for nodes */
/**  FUNCTION:  HMM_PROFILE_Reuse()
 *   SYNOPSIS:
 */
void HMM_PROFILE_SetModel_Length(HMM_PROFILE* prof,
                                 int length) {
  /* realloc memory if allocated length is less than new length */
  if (prof->Nalloc < length) {
    prof->hmm_model = ERROR_realloc(prof->hmm_model, (length + 1) * sizeof(HMM_NODE));
    prof->hmm_model_full = NULL;
    prof->Nalloc = length;
  }

  prof->N = length;
}

/* Set alphabet (DNA or AMINO ACID) for HMM_PROFILE */
/**  FUNCTION:  HMM_PROFILE_Reuse()
 *   SYNOPSIS:
 */
void HMM_PROFILE_SetAlphabet(HMM_PROFILE* prof,
                             char* alph_name) {
  prof->alph = STR_Create(alph_name);
  prof->alph = STR_ToLower(prof->alph);

  if (STR_Compare(alph_name, "amino") == 0) {
    prof->alph_leng = NUM_AMINO;
  }
  elif (STR_Compare(alph_name, "dna") == 0) {
    prof->alph_leng = NUM_DNA;
  }
  else {
    fprintf(stderr, "ERROR: Invalid alphabet type: %s\n", alph_name);
    ERRORCHECK_exit(EXIT_FAILURE);
  }
}

/* Determine the consensus sequence (highest likelihood output) using HMM_PROFILE */
/**  FUNCTION:  HMM_PROFILE_Reuse()
 *   SYNOPSIS:
 */
void HMM_PROFILE_SetConsensus(HMM_PROFILE* prof) {
  float best_val;
  float new_val;
  char best_amino;
  HMM_NODE curr_node;

  /* clear pre-existing consensus */
  VECTOR_CHAR_Reuse(prof->consensus);
  VECTOR_CHAR_GrowTo(prof->consensus, prof->N + 1);
  VECTOR_CHAR_SetSize(prof->consensus, prof->N + 1);

  /* find consensus sequence */
  for (int i = 1; i <= prof->N; i++) {
    best_val = -INF;
    curr_node = prof->hmm_model[i];
    /* for each residue in the consensus sequence, choose the most likely to be emitted */
    for (int j = 0; j < NUM_AMINO; j++) {
      new_val = curr_node.match[j];
      if (best_val < new_val) {
        best_val = new_val;
        best_amino = AA[j];
      }
    }
    VECTOR_CHAR_Set(prof->consensus, i - 1, best_amino);
  }
  VECTOR_CHAR_Set(prof->consensus, prof->N, NULL_CHAR);
}

/**  FUNCTION:  HMM_PROFILE_GetConsensus()
 *   SYNOPSIS:  Return <consensus> sequence as string.
 *              Will build it now if necessary.
 */
STR HMM_PROFILE_GetConsensus(HMM_PROFILE* prof) {
  STR consensus;
  /* build consensus now if not already built */
  if (prof->is_consensus == false) {
    HMM_PROFILE_SetConsensus(prof);
  }
  consensus = VECTOR_CHAR_GetArray(prof->consensus);

  /* If string is empty, return NULL */
  if (STR_GetLength(consensus) <= 0) {
    return NULL;
  }
  return consensus;
}

/* Set Distribution Parameters for HMM_PROFILE */
/**  FUNCTION:  HMM_PROFILE_Reuse()
 *   SYNOPSIS:
 */
void HMM_PROFILE_SetDistribution_Params(HMM_PROFILE* prof,
                                        float param1,
                                        float param2,
                                        char* dist_name) {
  float* parPtr1 = NULL;
  float* parPtr2 = NULL;

  if (STR_Compare(dist_name, "MSV") == 0) {
    prof->msv_dist.param1 = param1;
    prof->msv_dist.param2 = param2;
  } else if (STR_Compare(dist_name, "VITERBI") == 0) {
    prof->viterbi_dist.param1 = param1;
    prof->viterbi_dist.param2 = param2;
  } else if (STR_Compare(dist_name, "FORWARD") == 0) {
    prof->forward_dist.param1 = param1;
    prof->forward_dist.param2 = param2;
  } else {
    fprintf(stderr, "Invalid distribution type: %s\n", dist_name);
    ERRORCHECK_exit(EXIT_FAILURE);
  }

  parPtr1 = &param1;
  parPtr2 = &param2;
}

/* Constrain model to range [tbeg, tend] */
/**  FUNCTION:  HMM_PROFILE_Reuse()
 *   SYNOPSIS:
 */
void HMM_PROFILE_SetSubmodel(HMM_PROFILE* prof, int t_beg, int t_end) {
  /* check if range is acceptable */
  int T = prof->N;
  if (t_beg < 0 && t_end < 0 && t_beg > T && t_end > T) {
    printf("Subset range (%d,%d) is outside range (%d,%d).\n", t_beg, t_end, 0, T);
    ERRORCHECK_exit(EXIT_FAILURE);
  }

  /* save full model temporarily */
  prof->hmm_model_full = prof->hmm_model;
  prof->N_full = prof->N;

  /* current submodel */
  prof->hmm_model = &(prof->hmm_model[t_beg]);
  prof->N = t_end - t_beg;
}

/* Unconstrain model to cover entire */
/**  FUNCTION:  HMM_PROFILE_Reuse()
 *   SYNOPSIS:
 */
void HMM_PROFILE_UnsetSubmodel(HMM_PROFILE* prof) {
  /* override temporary model with full model */
  prof->hmm_model = prof->hmm_model_full;
  prof->N = prof->N_full;

  /* clear temporary save */
  prof->hmm_model_full = NULL;
  prof->N_full = -1;
}

/* Output HMM_PROFILE to FILE POINTER */
/**  FUNCTION:  HMM_PROFILE_Reuse()
 *   SYNOPSIS:
 */
void HMM_PROFILE_Dump(HMM_PROFILE* prof,
                      FILE* fp) {
  int pad_1, pad_2;

  /* test for bad file pointer */
  if (fp == NULL) {
    const char* obj_name = "HMM_PROFILE";
    fprintf(stderr, "ERROR: Bad FILE POINTER for printing %s.\n", obj_name);
    ERRORCHECK_exit(EXIT_FAILURE);
    return;
  }

  int i, j;
  if (prof->consensus == NULL) {
    HMM_PROFILE_SetConsensus(prof);
  }

  pad_1 = 15;
  fprintf(fp, "\n");
  fprintf(fp, "=== HMM PROFILE ====================================\n");
  fprintf(fp, "%*s %s\n", pad_1, "NAME:", prof->name);
  fprintf(fp, "%*s %d\n", pad_1, "LENGTH:", prof->N);
  fprintf(fp, "%*s %d\n", pad_1, "ALLOC:", prof->Nalloc);
  fprintf(fp, "%*s %s\n", pad_1, "CONSENSUS:", prof->consensus->data);

  /* background model */
  pad_1 = 13;
  pad_2 = 4;
  fprintf(fp, "#%*s:%*s", pad_1, "FREQ", pad_1, "");
  for (int j = 0; j < NUM_AMINO; j++)
    fprintf(fp, "%9.4f ", prof->bg_model->freq[j]);
  fprintf(fp, "\n");

  fprintf(fp, "#%*s:%*s", pad_1, "COMPO", pad_1, "");
  for (int j = 0; j < NUM_AMINO; j++)
    fprintf(fp, "%9.4f ", prof->bg_model->compo[j]);
  fprintf(fp, "\n");

  fprintf(fp, "#%*s:%*s", pad_1, "INSERT", pad_1, "");
  for (int j = 0; j < NUM_AMINO; j++)
    fprintf(fp, "%9.4f ", prof->bg_model->insert[j]);
  fprintf(fp, "\n");

  fprintf(fp, "#%*s:%*s", pad_1, "TRANS", pad_1, "");
  for (int j = 0; j < NUM_TRANS_STATES; j++)
    fprintf(fp, "%9.4f ", prof->bg_model->trans[j]);
  fprintf(fp, "\n\n");

  /* position-specific probabilities */
  fprintf(fp, " %5s ", "");
  for (int j = 0; j < NUM_AMINO_PLUS_SPEC; j++) {
    fprintf(fp, "%9d ", j);
  }
  fprintf(fp, "\n");
  for (int i = 0; i < prof->N + 1; i++) {
    /* line 1: match emissions */
    fprintf(fp, " %5d ", i);
    for (int j = 0; j < NUM_AMINO_PLUS_SPEC; j++)
      fprintf(fp, "%9.4f ", prof->hmm_model[i].match[j]);
    fprintf(fp, "\n");

    /* line 2: insert emissions */
    fprintf(fp, "%7s", "");
    for (int j = 0; j < NUM_AMINO_PLUS_SPEC; j++)
      fprintf(fp, "%9.4f ", prof->hmm_model[i].insert[j]);
    fprintf(fp, "\n");

    /* line 3: transition probabilities */
    fprintf(fp, "%7s", "");
    for (int j = 0; j < NUM_TRANS_STATES; j++)
      fprintf(fp, "%9.4f ", prof->hmm_model[i].trans[j]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "\n");

  /* background frequencies */
  pad_1 = 13;
  fprintf(fp, "#%*s:\n", pad_1, "BACKGROUND");

  fprintf(fp, "#%*s:%7s", pad_1, "LOG", "");
  for (int i = 0; i < NUM_AMINO; i++)
    fprintf(fp, "%9.4f ", BG_MODEL_log[i]);
  fprintf(fp, "\n");

  fprintf(fp, "#%*s:%7s", pad_1, "ACTUAL", "");
  for (int i = 0; i < NUM_AMINO; i++)
    fprintf(fp, "%9.4f ", BG_MODEL[i]);
  fprintf(fp, "\n");

  fprintf(fp, "#%15s:\n", "SPECIAL");
  fprintf(fp, "%16s:\t%9.4f %9.4f\n",
          "E", prof->bg_model->spec[SP_E][SP_LOOP], prof->bg_model->spec[SP_E][SP_MOVE]);
  fprintf(fp, "%16s:\t%9.4f %9.4f\n",
          "N", prof->bg_model->spec[SP_N][SP_LOOP], prof->bg_model->spec[SP_N][SP_MOVE]);
  fprintf(fp, "%16s:\t%9.4f %9.4f\n",
          "C", prof->bg_model->spec[SP_C][SP_LOOP], prof->bg_model->spec[SP_C][SP_MOVE]);
  fprintf(fp, "%16s:\t%9.4f %9.4f\n",
          "J", prof->bg_model->spec[SP_J][SP_LOOP], prof->bg_model->spec[SP_J][SP_MOVE]);

  fprintf(fp, "//\n");
}

void HMM_FILE_Dump(HMM_PROFILE* prof,
                   FILE* fp) {
  /* hardcoded vars */
  char* format_name = "HMMER3/f [3.2.1 | June 2018]";
  char* alph = "amino";
  int nseq = 1;

  int head_pad = 5;
  int stat_pad = 21;
  int hmm_pad = 7;
  int cons_pad = 6;

  /* Header */
  fprintf(fp, "%s\n", format_name);
  fprintf(fp, "%*s %s\n", head_pad, "NAME", prof->name);
  fprintf(fp, "%*s %d\n", head_pad, "LENG", prof->N);
  fprintf(fp, "%*s %s\n", head_pad, "ALPH", alph);
  /* RF - y/n */
  /* MM - y/n */
  /* CONS - y/n */
  /* CS - y/n */
  /* MAP - y/n */
  fprintf(fp, "%*s %d\n", head_pad, "NSEQ", nseq);
  /* CKSUM - checksum */
  fprintf(fp, "%*s %8.4f %8.5f\n",
          stat_pad, "STAT LOCAL MSV", prof->msv_dist.param1, prof->msv_dist.param2);
  fprintf(fp, "%*s %8.4f %8.5f\n",
          stat_pad, "STAT LOCAL VITERBI", prof->viterbi_dist.param1, prof->viterbi_dist.param2);
  fprintf(fp, "%*s %8.4f %8.5f\n",
          stat_pad, "STAT LOCAL FORWARD", prof->forward_dist.param1, prof->forward_dist.param2);

  /* hmm */
  fprintf(fp, "HMM          A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y\n");
  fprintf(fp, "            m->m     m->i     m->d     i->m     i->i     d->m     d->d\n");
  /* background composition */
  /* emission probabilities */
  fprintf(fp, "%*s   ", -hmm_pad, "COMPO");
  for (int j = 0; j < NUM_AMINO; j++) {
    fprintf(fp, "%8.5f ", prof->bg_model->compo[j]);
  }
  fprintf(fp, "\n");
  /* background insert probabilities */
  fprintf(fp, "%*s   ", -hmm_pad, "");
  for (int j = 0; j < NUM_AMINO; j++) {
    fprintf(fp, "%8.5f ", prof->bg_model->insert[j]);
  }
  fprintf(fp, "\n");
  /* back transition proabilities */
  fprintf(fp, "%*s   ", -hmm_pad, "");
  for (int j = 0; j < NUM_TRANS_STATES; j++) {
    if (prof->bg_model->trans[j] == -INF || prof->bg_model->trans[j] == INF) {
      fprintf(fp, "%8s ", "*");
    } else {
      fprintf(fp, "%8.5f ", prof->bg_model->trans[j]);
    }
  }
  fprintf(fp, "\n");
  /* position-specific probabilities */
  for (int i = 1; i < prof->N + 1; i++) {
    /* emission probabilities */
    fprintf(fp, "%*d   ", -hmm_pad, i);
    for (int j = 0; j < NUM_AMINO; j++) {
      fprintf(fp, "%8.5f ", prof->bg_model->compo[j]);
    }
    /* consensus */
    // fprintf(stderr, "%*d %c %c %c %c\n", -cons_pad, *prof->consensus[i], "-", "-", "-" );
    fprintf(fp, "\n");
    /* background insert probabilities */
    fprintf(fp, "%*s   ", -hmm_pad, "");
    for (int j = 0; j < NUM_AMINO; j++) {
      fprintf(fp, "%8.5f ", prof->bg_model->insert[j]);
    }
    fprintf(fp, "\n");
    /* back transition proabilities */
    fprintf(fp, "%*s   ", -hmm_pad, "");
    for (int j = 0; j < NUM_TRANS_STATES; j++) {
      if (prof->bg_model->trans[j] == -INF || prof->bg_model->trans[j] == INF) {
        fprintf(fp, "%8s ", "*");
      } else {
        fprintf(fp, "%8.5f ", prof->bg_model->trans[j]);
      }
    }
    fprintf(fp, "\n");
  }
}
