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
const int num_pipelines = 8;
PIPELINE PIPELINES[] = {
   /* mmore-seqs pipelines */
   { "search",       mmore_pipeline,            5,    NULL },
   { "mmore-main",   mmore_main_pipeline,       3,    NULL },
   { "prep",         prep_pipeline,             3,    NULL },
   { "mmore-prep",   mmore_prep_pipeline,       1,    NULL },
   { "easy-search",  mmore_easysearch_pipeline, 3,    NULL },
   /* helper pipelines */
   { "index",        index_pipeline,            2,    NULL },  
   { "hmmbuild",     hmmbuild_pipeline,         1,    NULL },
   { "interactive",  interactive_pipeline,      0,    NULL },
   /* alternative search pipelines */
   { "null",         null_pipeline,             0,    NULL },
   { "generic",      generic_pipeline,          2,    NULL },
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
char STATE_CHAR[] = {
   'M',
   'I',
   'D',
   'E',
   'N',
   'J',
   'C',
   'B',
   'S',
   'T',
   'X', 
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

/* Map file extensions to FILETYPE */
STR_TO_INT FILETYPE_EXTS[] = {
   { ".hmm",      FILE_HMM },     
   { ".fasta",    FILE_FASTA },
   { ".fa",       FILE_FASTA },
   { ".msa",      FILE_MSA },
   { ".hhm",      FILE_HHM },
   { ".mm_msa",   FILE_MM_MSA },
   { ".smmdb",    FILE_MMDB_S },
   { ".pmmdb",    FILE_MMDB_P }
};
const int NUM_FILETYPE_EXTS = 8;

int 
FILETYPE_EXT_Get( STR   filetype_ext )
{
   for (int i = 0; i < NUM_FILETYPE_EXTS; i++) {
      STR_TO_INT filetype_map = FILETYPE_EXTS[i];
      if ( STR_Compare( filetype_ext, filetype_map.s ) == 0 ) {
         return filetype_map.i;
      }
   }
}

/* Map FILETYPES to keyword names */
STR_TO_INT FILETYPE_NAMES[] = {
   { "NULL",      FILE_NULL },
   { "HMM",       FILE_HMM },
   { "FASTA",     FILE_FASTA },
   { "MSA",       FILE_MSA },
   { "MM_MSA",    FILE_MM_MSA },
   { "HHM",       FILE_HHM },
   { "MMDB",      FILE_MMDB },
   { "SMMDB",     FILE_MMDB_S },
   { "PMMDB",     FILE_MMDB_P }
};
const int NUM_FILETYPE_NAMES = 9;

STR 
FILETYPE_NAME_Get( FILETYPE   filetype_id )
{
   for (int i = 0; i < NUM_FILETYPE_NAMES; i++) {
      STR_TO_INT filetype_map = FILETYPE_NAMES[i];
      if ( filetype_id == filetype_map.i ) {
         return filetype_map.s;
      }
   }
   return NULL;
}

/* debugging data structures */
DEBUG_KIT*  debugger;

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