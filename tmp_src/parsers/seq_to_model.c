/*******************************************************************************
 *  FILE:      seq_to_model.c
 *  PURPOSE:   Converts a SEQUENCE to HMM_PROFILE.
 *             Uses the Easel and HMMER libraries.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/* import easel */
#include "../easel/easel.h"
#include "../easel/esl_dirichlet.h"
#include "../easel/esl_hmm.h"
#include "../easel/esl_sq.h"
#include "../easel/esl_sqio.h"
#include "../easel/esl_scorematrix.h"
#include "../easel/esl_getopts.h"
#include "../easel/esl_vectorops.h"

/* local imports */
#include "objects/structs.h"
#include "hmmer/p7_structs.h"
#include "hmmer/p7_funcs.h"
#include "objects/hmm_profile.h"
#include "objects/sequence.h"
#include "utilities/utility.h"

/* header */
#include "seq_to_profile.h"

/* converts sequence to HMM profile model */
void SEQUENCE_to_HMM_PROFILE( SEQUENCE*      seq, 
                  HMM_PROFILE*   prof )
{
   /* extract info */
   char*          filename = seq->filename;

   /* conversion data structs */
   ESL_SQFILE*    fp  = NULL;    /* Easel's file pointer */
   ESL_ALPHABET*  abc = NULL;    /* Easel's alphabet (set to amino acid) */
   P7_BG*         bg  = NULL;    /* HMMER's background null model */
   ESL_SQ*        qsq = NULL;    /* Easel's sequence datatype */
   P7_HMM*        hmm = NULL;    /* HMMER's hmm datatype */
   P7_PROFILE*    gm  = NULL;    /* HMMER's generic unoptimized hmm profile model */

   /* initialize alphabet and null model */
   // abc   = esl_alphabet_Create(eslAMINO);
   // bg    = p7_bg_Create(abc);
   // bld   = p7_builder_Create(NULL, abc);

   /* open file pointer */
   // esl_sqfile_OpenDigital(abc, seq->filename, qformat, NULL, &fp);

   /* convert SEQUENCE to ESL_SQ (digital sequence) */
   // qsq      = esl_sq_CreateDigital(abc);
   // qstatus  = esl_sqio_Read(qfp, qsq)

   /* Build the model from sequence */
   // p7_SingleBuilder(bld, qsq, bg, NULL, NULL, &gm, NULL); /* bypass HMM - only need model */

   /* convert P7_HMM to HMM_PROFILE */
}