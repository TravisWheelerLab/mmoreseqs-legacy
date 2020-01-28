/*******************************************************************************
 *  @file structs.h
 *  @brief All Data Structures used by Cloud Search.
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/

#ifndef _STRUCTS_H
#define _STRUCTS_H

extern char AA[];
extern int AA_REV[];
extern char AA2[];
extern int AA2_REV[];
extern float BG_MODEL[];
extern float BG_MODEL_log[];

/* MATCH, INSERT, DELETE, SPECIAL DP MATRIX ACCESS MACROS */
#define ST_MX(mx,st,i,j)   (   mx[ (    st*(Q+1)*(T+1)) + ((i)*(T+1)) + (j) ])
#define MMX(i,j)           (st_MX[ (MAT_ST*(Q+1)*(T+1)) + ((i)*(T+1)) + (j) ])
#define IMX(i,j)           (st_MX[ (INS_ST*(Q+1)*(T+1)) + ((i)*(T+1)) + (j) ])
#define DMX(i,j)           (st_MX[ (DEL_ST*(Q+1)*(T+1)) + ((i)*(T+1)) + (j) ])
/* MATCH, INSERT, DELETE, SPECIAL ANTI-DIAG ACCESS MACROS (d = diag, i = offset) */
#define MMX3(i,j)          (st_MX3[ (MAT_ST*3*(T+1)) + ((i)*(T+1)) + (j) ])
#define IMX3(i,j)          (st_MX3[ (INS_ST*3*(T+1)) + ((i)*(T+1)) + (j) ])
#define DMX3(i,j)          (st_MX3[ (DEL_ST*3*(T+1)) + ((i)*(T+1)) + (j) ])


/* SPECIAL STATE MATRIX MACROS */
#define SP_MX(mx,sp,i)     (mx[ ((sp)*(Q+1)) + (i) ])
#define XMX(sp,i)          (sp_MX[ ((sp)*(Q+1)) + (i) ])

/* TEST MATRIX */
#define TMX(i,j)           (test_MX[ ((i)*(T+1)) + (j) ])

/* TRANSITION SCORE, SPECIAL TRANSITION SCORE, MATCH SCORE, INSERT SCORE MACROS */
#define TSC(j,tr)       (target->hmm_model[j].trans[tr])
// #define TSC(j,tr)       (target->bg_model->trans[tr])
#define XSC(sp,tr)      (target->bg_model->spec[sp][tr])
#define MSC(j,A)        (target->hmm_model[j].match[A])
#define ISC(j,A)        (target->hmm_model[j].insert[A])

/* DEBUG MACRO FOR RETREIVING VARIABLE NAME */
#define getName(var) #var

/* BASIC FUNCTION MACROS */
#define MAX(i,j)     (( (i) > (j) ? (i) : (j) ))
#define MIN(i,j)     (( (i) < (j) ? (i) : (j) ))
#define ABS(i)       (( (i) > (0) ? (i) : (-i) ))
/* check if two value are equal within tolerance */
#define CMP_TOL(i,j) (( fabs( (i) - (j) ) < tol ? 1 : 0 )) 

/* DEFINED CONSTANTS */
#define CONST_LOG2 0.69314718055994529
#define SCALE_FACTOR 1000
#define INF INFINITY

/* max ASCII value of x - 'A' for alphabet */
#define ALPHA_MAX 26
#define SUBMAT_SIZE ALPHA_MAX*ALPHA_MAX

/* Search modes. */
typedef enum {
   MODE_NONE        = 0,    /* NO APPLICATIONS */
   MODE_MULTILOCAL  = 1,    /* multihit local:  "fs" mode   */
   MODE_MULTIGLOCAL = 2,    /* multihit glocal: "ls" mode   */
   MODE_UNILOCAL    = 3,    /* unihit local: "sw" mode      */
   MODE_UNIGLOCAL   = 4     /* unihit glocal: "s" mode      */
} SEARCH_MODE;

#define Test_IsLocal(mode)  (mode == MODE_MULTILOCAL || mode == MODE_UNILOCAL)
#define Test_IsMulti(mode)  (mode == MODE_MULTILOCAL || mode == MODE_MULTIGLOCAL)

typedef enum {
   MODE_DIAG = 0,
   MODE_ROW  = 1
} EDG_MODE;

typedef enum {
   MODE_LINEAR = 0,
   MODE_QUAD   = 1,
   MODE_NAIVE  = 2
} FWDBCK_MODE;

typedef enum {
   AMINO, DNA
} ALPHABET;

typedef struct {
   int st;
   int end;
} RANGE;

typedef struct {
   float param1;
   float param2;
} DIST_PARAM;

typedef struct {
   int i; /* row index */
   int j; /* col index */
} COORDS;

/* for debugging */
// const char* STATE_NAMES[] = {
//    "M_ST",
//    "I_ST",
//    "D_ST",
//    "E_ST",
//    "N_ST",
//    "J_ST",
//    "C_ST",
//    "B_ST",
//    "S_ST",
//    "T_ST",
//    "X_ST",
// };

typedef enum {
   M_ST = 0,
   I_ST = 1,
   D_ST = 2,
   E_ST = 3, 
   N_ST = 4, 
   J_ST = 5, 
   C_ST = 6, 
   B_ST = 7, 
   S_ST = 8,
   T_ST = 9,
   X_ST = 10
} ALL_STATES;
#define NUM_ALL_STATES 9

typedef enum {
   MAT_ST = 0,
   INS_ST = 1,
   DEL_ST = 2
} NORMAL_STATES;
#define NUM_NORMAL_STATES 3

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
   xGC = 20, /* gap character */
   xNC = 21, /* non-residue character */
   xMC = 22  /* missing character */
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
   SP_E = 0, /*  */
   SP_N = 1, /*  */
   SP_J = 2, /*  */
   SP_C = 3, /*  */
   SP_B = 4, /*  */
} SPECIAL_STATES;
#define NUM_SPECIAL_STATES 5

typedef struct {
   float match[NUM_AMINO];/* match emission probabilities for each amino acid */
   float insert[NUM_AMINO]; /* insert emission probabilities for each amino acid */
   /* [0]A  [1]C  [2]D  [3]E  [4]F  [5]G  [6]H  [7]I  [8]K  [9]L
      [10]M [11]N [12]P [13]Q [14]R [15]S [16]T [17]V [18]W [19]Y */
   float trans[NUM_TRANS_STATES]; /* transition state probabilities (default same as COMPO) */
   /* [0]m->m [1]m->i [2]m->d [3]i->m [4]i->i [5]d->m [6]d->d [7]b->m */
} HMM_NODE;

typedef struct {
   float freq[NUM_AMINO];    /* hard-coded background residue frequencies for each amino acid */
   float compo[NUM_AMINO];   /* background residue frequencies of the given hmm model */
   float insert[NUM_AMINO];  /* insert emission probabilities for each amino acid (uniform across positions) */
   /* [0]A  [1]C  [2]D  [3]E  [4]F  [5]G  [6]H  [7]I  [8]K  [9]L
      [10]M [11]N [12]P [13]Q [14]R [15]S [16]T [17]V [18]W [19]Y */
   float trans[NUM_TRANS_STATES]; /* transition state probabilities (default same as COMPO) */
   /* [0]m->m [1]m->i [2]m->d [3]i->m [4]i->i [5]d->m [6]d->d [7]b->m */

   /* special state transitions probabilities */
   float spec[NUM_SPECIAL_STATES][NUM_SPECIAL_TRANS];
   /* [0]N [1]E [2]C [3]J */
   /* [0]LOOP  [1]MOVE */
   int num_J;
} HMM_BG;

typedef struct {
   char *filename;
   char *name, *acc, *desc; /* meta data */
   int leng; /* sequence length */

   char *alph; /* ordering of alphabet */
   int alph_leng; /* alphabet length: AMINO = 20, DNA = 4 */

   DIST_PARAM *msv_dist, *viterbi_dist, *forward_dist; /* distribution parameters for scoring */

   HMM_BG *bg_model; /* background composition */
   HMM_NODE *hmm_model; /* array of position specific probabilities */

   int isLocal, isMultihit; /* profile settings */
} HMM_PROFILE;

typedef struct {
   char *filename;
   HMM_PROFILE profs;
   int num_profs;
} HMM_DB;

typedef struct {
   char *filename;
   char *name;
   int leng;
   char *alph;
   char *seq;
} SEQ;

typedef struct {
   char *filename;
   SEQ *seqs;
   int num_seqs;
} SEQ_DB;

typedef struct {
   char *filename;
   float *scores;
} SUBMAT;

typedef struct {
   float alpha; 
   int beta;
   int search_mode;
   char *target_hmm_file, *query_fasta_file;
   RANGE *range_query, *range_target;
   char *outfile_name;
   // FILE *outfile;
} ARGS;

typedef struct {
   float viterbi_sc;
   float fwd_sc, bck_sc;
   float cloud_fwd_sc, cloud_bck_sc;
   float perc_cells, perc_window;  
} SCORES;

typedef struct {
   SEQ *query;
   HMM_PROFILE *target;
   RANGE *query_range, *target_range;
   float score;
} HIT;

typedef struct {
   char *filename;
   int N;        /* current length */
   int size;     /* allocated length */
   HIT *hits;
} RESULTS;

typedef struct {
   int i;
   int j;
   int st;
} TRACE;

typedef struct {
   int N;                     /* current length */
   int size;                  /* allocated length */
   int L,W;                   /* length, width dimensions of query/target window */
   COORDS first_m, last_m;    /* (i,j) position of the last match found */
   float sc_last_m;
   int start, end;            /* position in trace for first/last MID state */
   float sc_max;              /* max alignment score in Viterbi */
   TRACE *traces;
} TRACEBACK;

typedef struct {
   int diag;                  /* current anti-diagonal OR row */
   int lb;                    /* bottom-left (lower) edge-bound */
   int rb;                    /* top-right (upper) edge-bound */
} BOUND;

typedef struct {
   int N;                     /* current length */
   int size;                  /* allocated length */
   int mode;                  /* 
   COORDS start, end;         /* starting point */
   BOUND *bounds;               
} EDGEBOUNDS;


#endif /* _STRUCTS_H */
