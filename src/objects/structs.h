/*******************************************************************************
 *  FILE:      structs.h
 *  PURPOSE:   All Data Structures used by Cloud Search.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _STRUCTS_H
#define _STRUCTS_H

/* === MACRO FUNCTIONS === */

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

#define Test_IsLocal(mode)  (mode == MODE_MULTILOCAL || mode == MODE_UNILOCAL)
#define Test_IsMulti(mode)  (mode == MODE_MULTILOCAL || mode == MODE_MULTIGLOCAL)

/* === ENUMERATIONS === */

/* Pipeline Modes */
typedef enum {
   PIPELINE_NULL     = 0,
   PIPELINE_MAIN     = 1,
   PIPELINE_TEST     = 2,
   PIPELINE_TIME     = 3,
   PIPELINE_MMSEQS   = 4,
   PIPELINE_INDEX    = 5,
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

typedef struct {
   float    viterbi_sc;
   float    fwd_sc;
   float    bck_sc;
   float    cloud_fwd_naive_sc; 
   float    cloud_bck_naive_sc;
   float    cloud_fwd_quad_sc; 
   float    cloud_bck_quad_sc;
   float    cloud_fwd_sc; 
   float    cloud_bck_sc;
   /* statistics */
   float    perc_cells; 
   float    perc_window;  
} SCORES;

/* === OBJECT STRUCTS == */
/* NOTE: These objects have their own .c files. */

typedef struct {
   int*     data;             /* array of data type */
   int      N;                /* current length of array in use */
   int      Nalloc;           /* current length of array allocated */
}  VECTOR_INT;

typedef struct {
   int      id;               /* current anti-diagonal OR row */
   int      lb;               /* bottom-left (lower) edge-bound */
   int      rb;               /* top-right (upper) edge-bound */
} BOUND;

typedef struct {
   int            N;          /* current size of bounds array */
   int            Nalloc;     /* allocated size of bounds array */
   VECTOR_INT*    ids;        /* identity of each row/diag */
   VECTOR_INT*    heads;      /* indexes the heads of each row/diag */
   BOUND*         bounds;     /* list of bounded ranges along a row/diag */
} EDGEBOUNDS;

typedef struct {
   int         i;             /* index in query */
   int         j;             /* index in target */
   int         st;            /* state at index */
} TRACE;

typedef struct {
   int         N;             /* current length */
   int         Nalloc;        /* allocated length */
   int         beg;           /* position in trace for first MID state */
   int         end;           /* position in trace for last MID state */
   TRACE*      traces;        /* list of all (state,i,j) TRACES in ALIGNMENT */
} ALIGNMENT;

typedef struct {
   time_t   start;            /* */
   time_t   stop;             /* */
   time_t   duration;         /* */

   float    N;                /* */
   float    Nalloc;           /* */
   float*   stamps;           /* */
} CLOCK;

typedef enum {
   AMINO    = 0,              /* Protein Alphabet */
   DNA      = 1,              /* DNA {ACGT} Alphabet */
} ALPHABET;

typedef struct {
   float    param1;
   float    param2;
} DIST_PARAM;

typedef struct {
   /* POSITION-SPECIFIC NORMAL STATE PROBABILITIES */
   /* match emission probabilities for each amino acid */
   float    match[NUM_AMINO];
   /* insert emission probabilities for each amino acid */  
   float    insert[NUM_AMINO];
   /* transition state probabilities (default same as COMPO) */
   float    trans[NUM_TRANS_STATES];   
} HMM_NODE;

typedef struct {
   /* BACKGROUND PROBABILITIES */
   /* hard-coded background residue frequencies for each amino acid */
   float    freq[NUM_AMINO]; 
   /* background residue frequencies of the given hmm model */
   float    compo[NUM_AMINO];
   /* insert emission probabilities for each amino acid (uniform across positions) */
   float    insert[NUM_AMINO];
   /* transition state probabilities (default same as COMPO) */
   float    trans[NUM_TRANS_STATES];   

   /* SPECIAL STATE PROBABILITIES */
   /* move, loop for special transition states */
   float    spec[NUM_SPECIAL_STATES][NUM_SPECIAL_TRANS];
   /* jump value for configuring HMM */
   int      num_J;
} HMM_BG;

typedef struct {
   int            N;                /* profile length (number of nodes) */
   int            alph_leng;        /* alphabet length: AMINO = 20, DNA = 4 */
   /* profile settings */
   int            isLocal;          /* */
   int            isMultihit;       /* */  
   /* file location */
   char*          filepath;         /* path to the file containing hmm */
   int            offset;           /* offset within file to hmm */
   /* meta data listed in hmm file */
   char*          name;             /* */
   char*          acc;              /* */
   char*          desc;             /* */
   char*          alph;             /* PROTEIN or DNA (Only PROTEIN ACCEPTED ) */;      
   /* distribution parameters for scoring */
   DIST_PARAM*    msv_dist;         /* Parameters for the Distribution for Ungapped Viterbi Scores */
   DIST_PARAM*    viterbi_dist;     /* Parameters for the Distribution for Viterbi Scores */
   DIST_PARAM*    forward_dist;     /* Parameters for the Distribution for Fwd/Bck Scores */
   /* models */
   HMM_BG*        bg_model;         /* background composition */
   HMM_NODE*      hmm_model;        /* array of position specific probabilities */
} HMM_PROFILE;

typedef struct {
   BOUND*   data;    /* array of data type */
   int      N;       /* current length of array in use */
   int      Nalloc;  /* current length of array allocated */
}  VECTOR_BOUND;

/* VECTOR struct */
typedef struct {
   TRACE*   data;     /* array of data type */
   int      N;        /* current length of array in use */
   int      Nalloc;   /* current length of array allocated */
}  VECTOR_TRACE;

typedef struct {
   int      N;         /* length of sequence */
   char*    filename;  /* filename of sequence */
   char*    name;      /* */
   char*    alph;      /* */
   char*    seq;       /* */
} SEQUENCE;

typedef struct {
   int      R;       /* number of columns = length of query */ 
   int      C;       /* number of rows = number of special states */
   int      Nalloc;  /* flat length of matrix = rows x cols */
   float*   data;    /* */
} MATRIX_2D;

typedef struct {
   int      R;       /* number of rows = length of query */
   int      C;       /* number of columns = length of target  */
   int      N;       /* number of 3rd dim = number of normal states */
   int      Nalloc;  /* */
   float*   data;    /* */
} MATRIX_3D;

/* === GLOBAL STRUCTS === */

typedef struct {
   /* file paths */
   char*    target_filepath;        /* filepath to target (hmm) file */
   char*    query_filepath;         /* filepath to query (fasta) file */
   /* index paths */
   char*    target_indexpath;       /* index filepath for quick access of target (hmm) file */
   char*    query_indexpath;        /* index filepath for quick access of query (fasta) file */
   /* mmseqs path */
   char*    hits_filepath;          /* filepath to output of MMSEQS pipeline  */
   /* output path */
   char*    output_filepath;        /* filepath to output results to; "!stdout" => stdout */

   /* cloud search tuning vars */
   float    alpha;                  /* cloud search: x-drop pruning ratio */
   int      beta;                   /* cloud search: number of antidiag passes before pruning  */

   /* pipeline options */
   int      pipeline_mode;          /* which workflow pipeline to use */
   int      verbosity_mode;         /* levels of verbosity */
   int      search_mode;            /* alignment search mode */

   /* if viterbi is precomputed, gives starting and ending coords (single result) */
   COORDS   beg;                    /* beginning coordinates of viterbi alignment */
   COORDS   end;                    /* ending coordinates of viterbi alignment */

   /* target/query metadata */
   int      target_filetype;        /* enumerated FILETYPE of target file */
   int      query_filetype;         /* enumerated FILETYPE of query file */

   int      target_fileno;          /*  */
   int      query_fileno;           /*  */ 
} ARGS;


typedef struct {
   float    load_hmm;

   float    viterbi; 
   float    traceback;

   float    fwd; 
   float    bck;

   float    cloud_fwd; 
   float    cloud_bck;

   float    merge;
   float    reorient;

   float    bound_fwd;
   float    bound_bck;
} TIMES;


typedef struct {
   char*       name;          /* Name of HMM/FASTA in file */
   long        offset;        /* Positional offset of HMM/FASTA into file */
} F_INDEX_NODE;

typedef struct {
   int            N;          /* Number of location index nodes used */
   int            Nalloc;     /* Number of location index nodes allocated */ 
   F_INDEX_NODE*  nodes;      /* List of nodes for location of each HMM/FASTA in file */

   char*          indexpath;  /* File path of index file (NULL if built and not loaded) */
   char*          filepath;   /* File path of file being indexed */
   int            filetype;   /* Type of file being indexed */

   char*          delim;      /* one of more delimiter of header fields */
   int            name_field; /* index of header field containing name */ 
   int            isSorted;   /* Whether the index nodes list has been sorted */
} F_INDEX;


typedef struct {
   int      num_args;      /* number of arguments */
   char*    name;          /* name of flag o */
   char*    long_flag;     /* long "--" flag */
   char*    short_flag;    /* single character "-" flag */
   char*    desc;          /* description of flag */
} FLAG_CMD;


typedef struct {
   char*    target_name;
   char*    query_name;

   COORDS   beg;
   COORDS   end;

   float    perc_id;
   int      aln_len;
   int      mismatch;
   int      gap_openings;

   int      query_start;
   int      query_end;
   int      target_start;
   int      target_end;

   double   e_value;
   int      bit_score;
} RESULT;

typedef struct {
   int      N;
   int      Nalloc;
   RESULT*  data;
   char*    filepath;
} RESULTS;


/* === GLOBAL VARIABLES === */

extern char*   STATE_NAMES[];
extern char*   STATE_FULL_NAMES[];

extern char*   PIPELINE_NAMES[];
extern void    (*PIPELINES[])(ARGS*);

extern char*   MODE_NAMES[];
extern char*   VERBOSITY_NAMES[];

extern char*   ALPHABET_NAMES[];
extern int     ALPHABET_LENGTHS[];

extern char*   FILE_TYPE_NAMES[];
extern int     FILE_TYPE_MAP[];

/* alphabetically-ordered amino lookup and reverse lookup tables */
extern char    AA[];
extern int     AA_REV[];
/* substitution matrix-ordered amino lookup and reverse lookup tables */
extern char    AA2[];
extern int     AA2_REV[];
/* background frequencies of null model, normal and log space */
extern float   BG_MODEL[];
extern float   BG_MODEL_log[];

#endif /* _STRUCTS_H */
