/*******************************************************************************
 *  FILE:      structs.h
 *  PURPOSE:   Enumerated Types and Constants used by Cloud Search.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _STRUCTS_ENUMS_H
#define _STRUCTS_ENUMS_H

/* === VERSION === */
#define BUILD_VERSION   0.1
#define BUILD_NAME      "alpha"
#define BUILD_DATE      __DATE__

/* === CONSTANTS === */
#define CONST_LOG2      0.69314718055994529     /* natural log: ln(2) */
#define SCALE_FACTOR    1000                    /* scaling factor for summing logrithm */
#define INF             INFINITY                /*  */
#define INT_MIN         -2147483648             /* min value of integer */

/* === OUTPUT PIPES === */
#define STDOUT          "/dev/stdout"
#define STDERR          "/dev/stderr"
#define DEVNULL         "/dev/null"
#define DEBUGOUT        "DEBUG.log"
#define DEBUG_VIZ       "DEBUG.viz"

#define debugout        debugger->dbfp
 
/* === ENUMERATIONS === */

/* commandline flags */
#define NUM_FLAG_CMDS 11

/* Method of sort for a list */
typedef enum{
   SORT_NONE,
   SORT_ID,
   SORT_NAME,
} SORT_TYPE;

/* Flag for which models in file are to be loaded */
typedef enum {
   LOAD_NONE,
   LOAD_BY_ID,       /* load by id (counting from start of file) */
   LOAD_BY_OFFSET,   /* laod by offset (referenced from index file) */
   LOAD_FIRST        /* load the first model in file */
} LOAD_TYPE;

/* Status Flags (for function returns) */
typedef enum {
   STATUS_SUCCESS,
   STATUS_FAILURE
} STATUS_FLAGS;
#define NUM_STATUS_FLAGS 2

/* Error Flags */
typedef enum {
   ERROR_NONE,
   ERROR_UNKNOWN,
   ERROR_MALLOC,
   ERROR_REALLOC,
   ERROR_FILE_IO,
   ERROR_UNSUPPORTED_FUNC,
} ERROR_FLAGS;
#define NUM_ERROR_FLAGS 6

/* Number Format of HMM_PROFILE */
typedef enum {
   PROF_FORMAT_REAL      = 0,
   PROF_FORMAT_NEGLOG    = 1,
   PROF_FORMAT_LOGODDS   = 2
} PROF_FORMAT;
#define NUM_PROF_FORMAT 3

/* Data Type */
typedef enum {
   DATATYPE_NONE     = 0,
   DATATYPE_INT      = 1,
   DATATYPE_FLOAT    = 2,
   DATATYPE_STRING   = 3,
   DATATYPE_BOOL     = 4,
} DATATYPES;
#define NUM_DATATYPES 5

/* Pipeline Modes */
typedef enum {
   PIPELINE_NULL     = 0,
   PIPELINE_MAIN     = 1,
   PIPELINE_TEST     = 2,
   PIPELINE_TIME     = 3,
   PIPELINE_MMSEQS   = 4,
   PIPELINE_INDEX    = 5,
   PIPELINE_UTEST    = 6,
   PIPELINE_VIZ      = 7,
} PIPELINE_MODE;
#define NUM_PIPELINE_MODES 8

/* Verbosity Modes */
typedef enum {
   VERBOSE_NONE   = 0,
   VERBOSE_LOW    = 1,
   VERBOSE_HIGH   = 2,
   VERBOSE_ALL    = 3
} VERBOSE_MODE;
#define NUM_VERBOSITY_MODES 3

/* select which targets and which queries to perform search against */
typedef enum {
   SELECT_NONE,            /* NO SEARCHES */
   SELECT_ALL_V_ALL,       /* Search all targets vs all queries */
   SELECT_FIRST_V_FIRST,   /* Search first target vs first query in file */
   SELECT_MMSEQS_LIST,     /* Search list of mmseqs hitlist (for mmseqs+) */
   SELECT_NAME_LIST,       /* Search list of names from targets/queries */
   SELECT_ID_LIST          /* Search list of ids from targets/queries */
} SELECT_SEARCH;

/* Search modes (cloud search only supports uniglocal) */
typedef enum {
   MODE_NULL        = 0,    /* NO APPLICATIONS */
   MODE_MULTILOCAL  = 1,    /* multihit local:  "fs" mode   */
   MODE_MULTIGLOCAL = 2,    /* multihit glocal: "ls" mode   */
   MODE_UNILOCAL    = 3,    /* unihit local: "sw" mode      */
   MODE_UNIGLOCAL   = 4,    /* unihit glocal: "s" mode      */
} SEARCH_MODE;

/* Flags whether EDGEBOUNDS are stored row-wise or antidiagonal-wise */
typedef enum {
   EDG_NONE,
   EDG_DIAG,
   EDG_ROW,
} EDG_MODE;

/* Flags whether using linear, quadratic, or naive algorithm for cloud search (TESTING) */
typedef enum {
   ALG_LINEAR = 0,
   ALG_QUAD   = 1,
   ALG_NAIVE  = 2
} ALG_MODE;

/* Flags whether output file should be overwritten or appended to */
typedef enum {
   MODE_OVERWRITE = 0,
   MODE_APPEND    = 1,
} WRITE_MODE;

/* Flags the filetype of incoming file */
typedef enum {
   FILE_NULL   = 0,
   FILE_HMM    = 1,
   FILE_FASTA  = 2,
} FILE_TYPE;
#define NUM_FILE_TYPES 3
#define NUM_FILE_EXTS  3

/*  */
typedef enum {
   M_ST = 0,   /* MATCH STATE */
   I_ST = 1,   /* INSERT STATE */
   D_ST = 2,   /* DELETE STATE */
   E_ST = 3,   /* END STATE */
   N_ST = 4,   /* NEW STATE */
   J_ST = 5,   /* JUMP STATE */
   C_ST = 6,   /* TERMINAL STATE */
   B_ST = 7,   /* BEGIN STATE */
   S_ST = 8,   /*  */
   T_ST = 9,   /*  */
   X_ST = 10   /*  */
} ALL_STATES;
#define NUM_ALL_STATES 9

/* Normal States that map into dynamic programming matrix */
typedef enum {
   MAT_ST = 0,    /* MATCH STATE */
   INS_ST = 1,    /* INSERT STATE */
   DEL_ST = 2     /* DELETE STATE */
} NORMAL_STATES;
#define NUM_NORMAL_STATES 3

/* Special States that map into dynamic programming matrix */
typedef enum {
   SP_E = 0,   /* END STATE */
   SP_N = 1,   /* NEW STATE */
   SP_J = 2,   /* JUMP STATE */
   SP_C = 3,   /* TERMINAL STATE */
   SP_B = 4,   /* BEGIN STATE */
} SPECIAL_STATES;
#define NUM_SPECIAL_STATES 5

/* Transition States that map into HMM_PROFILE */
typedef enum {
   M2M = 0,
   M2I = 1,
   M2D = 2,
   I2M = 3,
   I2I = 4,
   D2M = 5,
   D2D = 6,
   B2M = 7
} TRANS_STATES;
#define NUM_TRANS_STATES 8

typedef enum {
   AMINO_A = 0,
   AMINO_C = 1,
   AMINO_D = 2,
   AMINO_E = 3,
   AMINO_F = 4,
   AMINO_G = 5,
   AMINO_H = 6,
   AMINO_I = 7,
   AMINO_K = 8,
   AMINO_L = 9,
   AMINO_M = 10,
   AMINO_N = 11,
   AMINO_P = 12,
   AMINO_Q = 13,
   AMINO_R = 14,
   AMINO_S = 15,
   AMINO_T = 16,
   AMINO_V = 17,
   AMINO_W = 18,
   AMINO_Y = 19,
   /* non-amino special chars */
   xGC     = 20, /* gap character */
   xNC     = 21, /* non-residue character */
   xMC     = 22  /* missing character */
} AMINOS;
#define NUM_AMINO 20

typedef enum {
   DNA_A = 0,
   DNA_C = 1,
   DNA_G = 2,
   DNA_T = 3,
} DNAS;
#define NUM_DNA 4

typedef enum {
   SP_LOOP = 0,
   SP_MOVE = 1,
} SPECIAL_TRANS;
#define NUM_SPECIAL_TRANS 2

typedef enum {
   ALPH_NULL   = 0,     /* NULL Alphabet */
   ALPH_AMINO  = 1,     /* Protein Alphabet */
   ALPH_DNA    = 2,     /* DNA {ACGT} Alphabet */
} ALPH_TYPE;

#endif /* _STRUCTS_ENUMS_H */