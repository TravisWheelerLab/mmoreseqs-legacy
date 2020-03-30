/*******************************************************************************
 *  FILE:      structs.h
 *  PURPOSE:   Enumerated Types and Constants used by Cloud Search.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _STRUCTS_ENUMS_H
#define _STRUCTS_ENUMS_H

/* Data Type */
typedef enum {
   DATATYPE_NONE     = 0,
   DATATYPE_INT      = 1,
   DATATYPE_FLOAT    = 2,
   DATATYPE_STRING   = 3,
} DATATYPES;
#define NUM_DATATYPES 4

/* Pipeline Modes */
typedef enum {
   PIPELINE_NULL     = 0,
   PIPELINE_MAIN     = 1,
   PIPELINE_TEST     = 2,
   PIPELINE_TIME     = 3,
   PIPELINE_MMSEQS   = 4,
   PIPELINE_INDEX    = 5,
   PIPELINE_UTEST    = 6,
} PIPELINE_MODE;
#define NUM_PIPELINE_MODES 6

/* Verbosity Modes */
typedef enum {
   VERBOSE_NONE   = 0,
   VERBOSE_LOW    = 1,
   VERBOSE_HIGH   = 2,
} VERBOSITY_MODE;
#define NUM_VERBOSITY_MODES 3

/* Search modes */
typedef enum {
   MODE_NULL        = 0,    /* NO APPLICATIONS */
   MODE_MULTILOCAL  = 1,    /* multihit local:  "fs" mode   */
   MODE_MULTIGLOCAL = 2,    /* multihit glocal: "ls" mode   */
   MODE_UNILOCAL    = 3,    /* unihit local: "sw" mode      */
   MODE_UNIGLOCAL   = 4,    /* unihit glocal: "s" mode      */
} SEARCH_MODE;

typedef enum {
   MODE_DIAG = 0,
   MODE_ROW  = 1,
} EDG_MODE;

typedef enum {
   MODE_LINEAR = 0,
   MODE_QUAD   = 1,
   MODE_NAIVE  = 2
} FWDBCK_MODE;

typedef enum {
   MODE_OVERWRITE = 0,
   MODE_APPEND    = 1,
} WRITE_MODE;

typedef enum {
   FILE_NULL   = 0,
   FILE_HMM    = 1,
   FILE_FASTA  = 2,
} FILE_TYPE;
#define NUM_FILE_TYPES 3
#define NUM_FILE_EXTS  3

typedef struct {
   int beg;
   int end;
} RANGE;

typedef struct {
   int i; /* row index */
   int j; /* col index */
} COORDS;

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

typedef enum {
   MAT_ST = 0,    /* MATCH STATE */
   INS_ST = 1,    /* INSERT STATE */
   DEL_ST = 2     /* DELETE STATE */
} NORMAL_STATES;
#define NUM_NORMAL_STATES 3

typedef enum {
   SP_E = 0,   /* END STATE */
   SP_N = 1,   /* NEW STATE */
   SP_J = 2,   /* JUMP STATE */
   SP_C = 3,   /* TERMINAL STATE */
   SP_B = 4,   /* BEGIN STATE */
} SPECIAL_STATES;
#define NUM_SPECIAL_STATES 5

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

typedef struct {
   float vals[NUM_TRANS_STATES];
   /* [0]m->m  [1]m->i  [2]m->d  [3]i->m  [4]i->i  [5]d->m  [6]d->d */
} TRANS_PROB;

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
   AMINO    = 0,              /* Protein Alphabet */
   DNA      = 1,              /* DNA {ACGT} Alphabet */
} ALPHABET;


#endif /* _STRUCTS_ENUMS_H */