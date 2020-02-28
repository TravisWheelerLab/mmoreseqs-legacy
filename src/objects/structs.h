/*******************************************************************************
 *  FILE:      structs.h
 *  PURPOSE:   All Data Structures used by Cloud Search.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _STRUCTS_H
#define _STRUCTS_H

extern char* ALL_STATE_NAMES[];
extern char  AA[];
extern int   AA_REV[];
extern char  AA2[];
extern int   AA2_REV[];
extern float BG_MODEL[];
extern float BG_MODEL_log[];

/* MATCH, INSERT, DELETE, SPECIAL DP MATRIX ACCESS MACROS */
#define ST_MX(mx,st,i,j)   ( mx[ ( st*(Q+1)*(T+1) ) + ( (i)*(T+1) ) + (j) ] )
#define MMX(i,j)           ST_MX( st_MX, MAT_ST, i, j)
#define IMX(i,j)           ST_MX( st_MX, INS_ST, i, j)
#define DMX(i,j)           ST_MX( st_MX, DEL_ST, i, j)
/* MATCH, INSERT, DELETE, SPECIAL ANTI-DIAG ACCESS MACROS (d = diag, i = offset) */
#define ST_MX3(mx,st,i,j)  ( mx[ ( st*3*((T+1)+(Q+1)) ) + ( (i)*((T+1)+(Q+1)) ) + (j) ] )
#define MMX3(i,j)          ST_MX3( st_MX3, MAT_ST, i, j)
#define IMX3(i,j)          ST_MX3( st_MX3, INS_ST, i, j)
#define DMX3(i,j)          ST_MX3( st_MX3, DEL_ST, i, j)

/* SPECIAL STATE MATRIX MACROS */
#define SP_MX(mx,sp,i)     (mx[ ((sp)*(Q+1)) + (i) ])
#define XMX(sp,i)          (sp_MX[ ((sp)*(Q+1)) + (i) ])

/* TEST MATRIX */
#define TMX(i,j)           (test_MX[ ((i)*(T+1)) + (j) ])

/* TRANSITION SCORE, SPECIAL TRANSITION SCORE, MATCH SCORE, INSERT SCORE MACROS */
#define TSC(j,tr)          (target->hmm_model[j].trans[tr])
// #define TSC(j,tr)          (target->bg_model->trans[tr])
#define XSC(sp,tr)         (target->bg_model->spec[sp][tr])
#define MSC(j,A)           (target->hmm_model[j].match[A])
#define ISC(j,A)           (target->hmm_model[j].insert[A])

/* DEBUG MACRO FOR RETREIVING VARIABLE NAME */
#define getName(var) #var

/* BASIC FUNCTION MACROS */
#define MAX(x,y)     (((x) > (y)) ? (x) : (y))
#define MIN(x,y)     (((x) < (y)) ? (x) : (y))
#define ABS(i)       (( (i) > (0) ? (i) : (-i) ))
/* check if two value are equal within tolerance */
#define CMP_TOL(i,j) (( fabs( (i) - (j) ) < tol ? 1 : 0 )) 

/* DEFINED CONSTANTS */
#define CONST_LOG2 0.69314718055994529
#define SCALE_FACTOR 1000
#define INF INFINITY
#define INT_MIN -2147483648

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
   MODE_OVERWRITE = 0,
   MODE_APPEND    = 1,
} WRITE_MODE;

typedef struct {
   int beg;
   int end;
} RANGE;

typedef struct {
   int i; /* row index */
   int j; /* col index */
} COORDS;

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
   SP_E = 0, /* END STATE */
   SP_N = 1, /* NEW STATE */
   SP_J = 2, /* JUMP STATE */
   SP_C = 3, /* TERMINAL STATE */
   SP_B = 4, /* BEGIN STATE */
} SPECIAL_STATES;
#define NUM_SPECIAL_STATES 5

typedef struct {
   float          viterbi_sc;
   float          fwd_sc;
   float          bck_sc;
   float          cloud_fwd_naive_sc; 
   float          cloud_bck_naive_sc;
   float          cloud_fwd_quad_sc; 
   float          cloud_bck_quad_sc;
   float          cloud_fwd_sc; 
   float          cloud_bck_sc;
   /* statistics */
   float          perc_cells; 
   float          perc_window;  
} SCORES;

/* === OBJECT STRUCTS == */
/* NOTE: Move to header files of their objects? */

typedef struct {
   int *data;        /* array of data type */
   int N;            /* current length of array in use */
   int Nalloc;       /* current length of array allocated */
}  VECTOR_INT;

typedef struct {
   int      id;                    /* current anti-diagonal OR row */
   int      lb;                    /* bottom-left (lower) edge-bound */
   int      rb;                    /* top-right (upper) edge-bound */
} BOUND;

typedef struct {
   int            N;                 /*  */
   int            Nalloc;            /*  */
   VECTOR_INT     *ids;              /* identity of each row/diag */
   VECTOR_INT     *heads;            /* indexes the heads of each row/diag */
   BOUND          *bounds;           /* list of bounded ranges along a row/diag */
} EDGEBOUNDS;

typedef struct {
   int         i;            /* index in query */
   int         j;            /* index in target */
   int         st;           /* state at index */
} TRACE;

typedef struct {
   int         N;            /* current length */
   int         Nalloc;       /* allocated length */
   int         beg;          /* position in trace for first MID state */
   int         end;          /* position in trace for last MID state */
   TRACE       *traces;      /* list of all (state,i,j) TRACES in ALIGNMENT */
} ALIGNMENT;

typedef struct {
   time_t   start;
   time_t   stop;
   time_t   duration;

   float    N;
   float    Nalloc;
   float*   stamps;
} CLOCK;

typedef enum {
   AMINO, 
   DNA
} ALPHABET;

typedef struct {
   float param1;
   float param2;
} DIST_PARAM;

typedef struct {
   /* NORMAL STATE PROBABILITIES */
   float          match[NUM_AMINO];          /* match emission probabilities for each amino acid */
   float          insert[NUM_AMINO];         /* insert emission probabilities for each amino acid */
   /* [0]A  [1]C  [2]D  [3]E  [4]F  [5]G  [6]H  [7]I  [8]K  [9]L
      [10]M [11]N [12]P [13]Q [14]R [15]S [16]T [17]V [18]W [19]Y */
   float          trans[NUM_TRANS_STATES];   /* transition state probabilities (default same as COMPO) */
   /* [0]m->m [1]m->i [2]m->d [3]i->m [4]i->i [5]d->m [6]d->d [7]b->m */
} HMM_NODE;

typedef struct {
   float          freq[NUM_AMINO];           /* hard-coded background residue frequencies for each amino acid */
   float          compo[NUM_AMINO];          /* background residue frequencies of the given hmm model */
   float          insert[NUM_AMINO];         /* insert emission probabilities for each amino acid (uniform across positions) */
   /* [0]A  [1]C  [2]D  [3]E  [4]F  [5]G  [6]H  [7]I  [8]K  [9]L
      [10]M [11]N [12]P [13]Q [14]R [15]S [16]T [17]V [18]W [19]Y */
   float          trans[NUM_TRANS_STATES];   /* transition state probabilities (default same as COMPO) */
   /* [0]m->m [1]m->i [2]m->d [3]i->m [4]i->i [5]d->m [6]d->d [7]b->m */

   /* SPECIAL STATE PROBABILITIES */
   float          spec[NUM_SPECIAL_STATES][NUM_SPECIAL_TRANS];
   /* [0]N [1]E [2]C [3]J */
   /* [0]LOOP  [1]MOVE */
   int            num_J;
} HMM_BG;

typedef struct {
   int            N;                   /* profile length (number of nodes) */
   int            alph_leng;           /* alphabet length: AMINO = 20, DNA = 4 */
   /* profile settings */
   int            isLocal; 
   int            isMultihit;     
   /* */
   char           *filename;
   /* meta data */
   char           *name; 
   char           *acc; 
   char           *desc;
   char           *alph;               /* alphabet type () */;      
   /* distribution parameters for scoring */
   DIST_PARAM     *msv_dist; 
   DIST_PARAM     *viterbi_dist; 
   DIST_PARAM     *forward_dist; 
   /* models */
   HMM_BG         *bg_model;           /* background composition */
   HMM_NODE       *hmm_model;          /* array of position specific probabilities */
} HMM_PROFILE;

typedef struct {
   int R, C;
   int Nalloc;
   float *data;
} MATRIX;

typedef struct {
   BOUND *data;      /* array of data type */
   int N;            /* current length of array in use */
   int Nalloc;       /* current length of array allocated */
}  VECTOR_BOUND;

/* VECTOR struct */
typedef struct {
   TRACE       *data;        /* array of data type */
   int         N;            /* current length of array in use */
   int         Nalloc;       /* current length of array allocated */
}  VECTOR_TRACE;

typedef struct {
   int    N;
   char*  filename;
   char*  name;
   char*  alph;
   char*  seq;
} SEQUENCE;

typedef struct {
   int     R; 
   int     C;
   int     Nalloc;
   float*  data;
} MATRIX_2D;

typedef struct {
   int     R; 
   int     C;
   int     N;
   int     Nalloc;
   float*  data;
} MATRIX_3D;

#endif /* _STRUCTS_H */
