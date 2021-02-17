/*******************************************************************************
 *     FILE:   structs.c
 *  PURPOSE:   All Data Structures and Global Static Structs used by Cloud Search.
 *
 *  AUTHOR:    Dave Rich
 *     BUG:    
 *******************************************************************************/

/* imports */
#include <time.h>
#include <stdio.h>
#include <stdbool.h>

/* local imports */
#include "../utilities/_utilities.h"
#include "_objects.h"
#include "../parsers/_parsers.h"
#include "../pipelines/_pipelines.h"

/* header import */
#include "structs.h"

/* Full Protein Alphabet with all Degen characters */
char ALPH_AMINO_CHARS[] = "ACDEFGHIKLMNPQRSTVWY-BJZOUX~";

/* Proteins in Alphabetical order */
char AA[] = { 
   'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
   'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
   /* special characters */
   '-', 
   /* degen characters */
   'X', 
   /* special characters */
   '*', '~'
};

/* Maps ASCII Letters to corresponding letters in AA[] */
int AA_REV[] = {
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1,  0, -1,  1,  2,  3, /* begin uppercase alphabet */
    4,  5,  6,  7, -1,  8, 9,  10, 11, -1,
   12, 13, 14, 15, 16, -1, 17, 18, 20, 19, /* X - unknown character */
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 
};

/* Null model background frequencies (normal space) (ported from HMMER) */
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

/* Null model background -log(x) frequencies (neg log space) (ported by HMMER) */
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

/* descriptors of all pipelines */
const int num_pipelines = 7;
PIPELINE PIPELINES[] = {
   { "null",         null_pipeline,          0 },
   { "generic",      generic_pipeline,       2 },
   { "mmore",        mmore_pipeline,         4 },
   { "mmore_main",   mmore_main_pipeline,    2 },
   { "index",        index_pipeline,         2 },
   { "hmmbuild",     hmmbuild_pipeline,      1 },
   { "interactive",  interactive_pipeline,   0 }
};

/* full names of the all states */
char* STATE_FULL_NAMES[] = { 
   "MATCH",
   "INSERT",
   "DELETE",
   "END",
   "INIT",
   "NEW",
   "JUMP",
   "INIT",
   "BEGIN",
   "START",
   "TERMINAL",
   "UNKNOWN",
};

/* Abbreviations of all states (for trace outputs) */
char* STATE_CHARS[] = {
   "M",
   "I",
   "D",
   "E",
   "N",
   "J",
   "C",
   "B",
   "S",
   "T",
   "X", 
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
   "High",
   "All",
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
   ".fa",
   ".msa",
   ".hhm",
   ".mm_msa"
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
   "MSA",
   "MM_MSA",
   "HHM"
};

/* command line flags and options */
int   num_flag_cmds = 11;
ARG_OPT COMMAND_OPTS[] = {
   /* name | num_args | data_type | arg_loc | arg_bool | long_flag | short_flag | desc */
   /* output formats */
   {  "OUTFILE",           1,    DATATYPE_INT,        NULL,    "--output",          "-o",    "Result output file destination [test_output/results.tsv]."  },
   {  "INFILES",           2,    DATATYPE_STRING,     NULL,    "--input",           "-i",    "Input files: {target,query} [test cases]."  },
   /* general options */
   {  "INDEX",             2,    DATATYPE_STRING,     NULL,    "--index",           "-x",    "Index files: {target,query} [builds on fly]."  },
   /* mmseqs options */
   {  "MMSEQS_TMP",        1,    DATATYPE_STRING,     NULL,    "--mmseqs-tmp",      NULL,    "MMseqs temp folder [null]."  },
   {  "MMSEQS_INPUT",      1,    DATATYPE_STRING,     NULL,    "--mmseqs-input",    NULL,    "MMseqs results file input [null]."  },
   {  "MMSEQS_LOOKUP",     2,    DATATYPE_STRING,     NULL,    "--mmseqs-lookup",   NULL,    "MMseqs lookup files: {target,query} [null]."  },
   /* mmore main options */
   {  "ALPHA",             1,    DATATYPE_FLOAT,      NULL,    "--alpha",           "-a",    "MMORE X-drop per antidiagonal pruning ratio [20.0]." },
   {  "BETA",              1,    DATATYPE_INT,        NULL,    "--beta",            "-b",    "MMore X-drop global " },
   {  "BETA",              1,    DATATYPE_INT,        NULL,    "--beta",            "-b",    "Number of passes of cloud search before pruning [5]." },
   
   {  "WINDOW",            4,    DATATYPE_INT,        NULL,    "--window",          "-w",    "Examine substring of query and target."  },
   {  "Q_RANGE",           2,    DATATYPE_INT,        NULL,    "--qrange",          NULL,    "Give range of ids in query file index to search [-1,-1]."  },
   {  "T_RANGE",           2,    DATATYPE_INT,        NULL,    "--trange",          NULL,    "Give range of ids in target file index to search [-1,-1]."  },
};

/* debugging data structures */
DEBUG_KIT*  debugger;

char* DATATYPE_NAMES[] = {
   "NONE",
   "INT",
   "FLOAT",
   "STR",
   "BOOL"
};


/* initializes single_builder scoring matrix used in HMM_PROFILE_From_Seq() */
/* SingleBuilder probability matrix (created by HMMER using p7_SingleBuilder, dependent on only bg frequencies) */
/* TODO: Needs to be freed */ 
SCORE_MATRIX* bld = NULL;


/* --- EXTERNAL EXECUTABLE/SCRIPT LOCATIONS --- */
char* ROOT_DIR                = MACRO_XSTR(PROJECT_LOC);
/* mmore-seqs script location */
char* MMSEQS_PLUS_SCRIPT      = MACRO_XSTR(PROJECT_LOC) "/scripts/workflows/mmseqs_plus.sh";
/* mmore-seqs script location */
char* MMSEQS_PLUS_EASY_SCRIPT = MACRO_XSTR(PROJECT_LOC) "/scripts/workflows/mmseqs_plus_easy.sh";
/* fasta-to-hmm converter script location */
char* FASTA_TO_HMM_SCRIPT     = MACRO_XSTR(PROJECT_LOC) "/scripts/workflows/convert_fasta_to_hmm.sh";
/* mmseqs binary location */
char* MMSEQS_BIN              = MACRO_XSTR(PROJECT_LOC) "/" MACRO_XSTR(MMSEQS_BIN_LOC);
/* hmmer 'hmmbuild' binary location */
char* HMMER_BIN               = MACRO_XSTR(PROJECT_LOC) "/" MACRO_XSTR(HMMER_BIN_LOC);
/* hmmer 'mmore' binary location */
char* MMORE_BIN               = MACRO_XSTR(PROJECT_LOC) "/" MACRO_XSTR(MMORE_BIN_LOC);

