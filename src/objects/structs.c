/*******************************************************************************
 *  @file structs.c
 *  @brief Data Structures for Cloud Search.
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

/* imports */
#include <time.h>
#include <stdio.h>
#include <stdbool.h>

/* header */
#include "structs.h"

/* local imports */
#include "utilities/utility.h"
#include "pipeline/pipeline.h"
#include "parsers/index_parser.h"

char ALPH_AMINO_CHARS[] = "ACDEFGHIKLMNPQRSTVWY-BJZOUX~";

/* Proteins in Alphabetical order */
char AA[] = { 
   'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
   'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' 
};

/* Maps ASCII Code to corresponding letters in AA2[] */
int AA_REV[] = {
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1,  0, -1,  1,  2,  3,
    4,  5,  6,  7, -1,  8, 9,  10, 11, -1,
   12, 13, 14, 15, 16, -1, 17, 18, -1, 19,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  
};

/* Null model background frequencies (normal space) used by HMMER */
double BG_MODEL[] = {
   0.0787945,     /* A */
   0.0151600,     /* C */
   0.0535222,     /* D */
   0.0668298,     /* E */
   0.0397062,     /* F */
   0.0695071,     /* G */
   0.0229198,     /* H */
   0.0590092,     /* I */
   0.0594422,     /* K */
   0.0963728,     /* L */
   0.0237718,     /* M */
   0.0414386,     /* N */
   0.0482904,     /* P */
   0.0395639,     /* Q */
   0.0540978,     /* R */
   0.0683364,     /* S */
   0.0540687,     /* T */
   0.0673417,     /* V */
   0.0114135,     /* W */
   0.0304133,     /* Y */
};

/* Null model background -log(x) frequencies (neg log space) used by HMMER */
double BG_MODEL_log[] = {
   2.540912081508539,     /* A */
   4.189094898767911,     /* C */
   2.927658757878447,     /* D */
   2.705606190131599,     /* E */
   3.226247932197823,     /* F */
   2.666326373355810,     /* G */
   3.775754113177160,     /* H */
   2.830061915029190,     /* I */
   2.822750867144509,     /* K */
   2.339531274855952,     /* L */
   3.739255274772411,     /* M */
   3.183542465385378,     /* N */
   3.030522495842527,     /* P */
   3.229838192658035,     /* Q */
   2.916961759390943,     /* R */
   2.683312711470047,     /* S */
   2.917499818784601,     /* T */
   2.697975620542613,     /* V */
   4.472958413679587,     /* W */
   3.492875266245181,     /* Y */
};

/* function pointers to all pipelines */
void (*PIPELINES[])(WORKER*) = {
   NULL,
   main_pipeline,
   test_pipeline,
   time_pipeline,
   mmseqs_pipeline,
   index_pipeline
};

/* pipeline names (for outputs) */
char* PIPELINE_NAMES[] = {
   "null",
   "main",
   "test",
   "time",
   "mmseqs",
   "index",
   "utest"
};

/* full names of the all states */
char* STATE_FULL_NAMES[] = { 
   "MATCH",
   "INSERT",
   "DELETE",
   "END",
   "",
   "JUMP",
   "BEGIN",
   "START",
   "",
   "UNKNOWN",
};

/* Abbreviations of all states (for trace outputs) */
char* STATE_NAMES[] = {
   "ST_M",
   "ST_I",
   "ST_D",
   "ST_E",
   "ST_N",
   "ST_J",
   "ST_C",
   "ST_B",
   "ST_S",
   "ST_T",
   "ST_X", 
};

/* Search Mode Names (for output) */
char* MODE_NAMES[] = { 
   "None", 
   "Multi-local", 
   "Multi-glocal", 
   "Uni-local", 
   "Uni-glocal",
}; 

/* Verbosity Level Names */
char* VERBOSITY_NAMES[] = {
   "None",
   "Low",
   "High"
};

/* Alphabet Names (for hmm files) */
char* ALPHABET_NAMES[] = {
   "AMINO",
   "DNA",
};

/* Alphabet Lengths */
int ALPHABET_LENGTHS[] = {
   20,      /* amino alphabet length */
   4        /* DNA alphabet length */
};

/* file extension types */
char* FILE_TYPE_EXTS[] = {
   ".hmm",     
   ".fasta",
   ".fa"
};

/* maps FILE_TYPE_NAMES to enum FILE_TYPE */
int FILE_TYPE_MAP[] = {
   FILE_HMM,
   FILE_FASTA,
   FILE_FASTA,
};

/* FILE_NAMES that map to FILE_TYPES */
char* FILE_TYPE_NAMES[] = {
   "NULL",
   "HMMER",
   "FASTA",
};

/* commandline arguments */
ARGS* args;

/* command line flags and options */
/* NOTE: update definition of NUM_FLAG_CMDS */
FLAG_CMD COMMAND_OPTS[] = {
   {  "ALPHA",    1,    DATATYPE_FLOAT,   NULL,    "--alpha",     "-a",    "X-drop pruning ratio." },
   {  "BETA",     1,    DATATYPE_INT,     NULL,    "--beta",      "-b",    "Number of passes of cloud search before pruning." },
   {  "WINDOW",   4,    DATATYPE_INT,     NULL,    "--window",    "-w",    "Examine substring of query and target."  },
   {  "Q_RANGE",  2,    DATATYPE_INT,     NULL,    "--qrange",    NULL,    "Give range of ids in query file index to search."  },
   {  "T_RANGE",  2,    DATATYPE_INT,     NULL,    "--trange",    NULL,    "Give range of ids in target file index to search."  },
   {  "SINGLE",   2,    DATATYPE_INT,     NULL,    "--single",    NULL,    "Give single {target,query} id pair to search." }
};

char* DATATYPE_NAMES[] = {
   "NONE",
   "INT",
   "FLOAT",
   "STR",
   "BOOL"
};
