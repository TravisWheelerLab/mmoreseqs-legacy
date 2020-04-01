/*******************************************************************************
 *  FILE:      structs.h
 *  PURPOSE:   All Data Structures used by Cloud Search.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _STRUCTS_H
#define _STRUCTS_H

/* === STDLIB DATA TYPES === */
#include <stdbool.h>

/* === MACRO FUNCTIONS === */
#include "objects/structs_macros.h"

/* === ENUMERATIONS === */
#include "objects/structs_enums.h"

/* === STRUCTS === */

/* */
typedef struct {
   int beg;
   int end;
} RANGE;

/* */
typedef struct {
   int i; /* row index */
   int j; /* col index */
} COORDS;

/* */
typedef struct {
   /* naive */
   float    cloud_fwd_naive_sc; 
   float    cloud_bck_naive_sc;
   /* quadratic */
   float    viterbi_quad_sc;
   float    fwd_quad_sc;
   float    bck_quad_sc;  
   float    cloud_fwd_quad_sc; 
   float    cloud_bck_quad_sc;
   /* linear */
   float    viterbi_sc;
   float    fwd_sc;
   float    bck_sc;
   float    cloud_fwd_sc; 
   float    cloud_bck_sc;
   /* statistics */
   float    perc_cells; 
   float    perc_window;  
} SCORES;

/* vector of integers */
typedef struct {
   int*     data;             /* array of data type */
   int      N;                /* current length of array in use */
   int      Nalloc;           /* current length of array allocated */
}  VECTOR_INT;

/* bound for single row or diagonal */
typedef struct {
   int      id;               /* current anti-diagonal OR row */
   int      lb;               /* bottom-left (lower) edge-bound */
   int      rb;               /* top-right (upper) edge-bound */
} BOUND;

/* set of bounds for cloud search space */
typedef struct {
   int            N;          /* current size of bounds array */
   int            Nalloc;     /* allocated size of bounds array */
   VECTOR_INT*    ids;        /* identity of each row/diag */
   VECTOR_INT*    heads;      /* indexes the heads of each row/diag */
   BOUND*         bounds;     /* list of bounded ranges along a row/diag */
} EDGEBOUNDS;

/* given cell of alignment */
typedef struct {
   int         i;             /* index in query */
   int         j;             /* index in target */
   int         st;            /* state at index */
} TRACE;

/* alignment for viterbi traceback */
typedef struct {
   int         N;             /* current length */
   int         Nalloc;        /* allocated length */
   int         beg;           /* position in trace for first MID state */
   int         end;           /* position in trace for last MID state */
   TRACE*      traces;        /* list of all (state,i,j) TRACES in ALIGNMENT */
} ALIGNMENT;

/* clock for timing events */
typedef struct {
   time_t   start;            /* */
   time_t   stop;             /* */
   time_t   duration;         /* */

   float    N;                /* */
   float    Nalloc;           /* */
   float*   stamps;           /* */
} CLOCK;

/* distribution parameters */ 
typedef struct {
   float    param1;
   float    param2;
} DIST_PARAM;

/* */
typedef struct {
   int      K;             /* size of unique alphabet */
   int      Kp;            /* size of total symbols: alphabet + special symbols  */
   char*    sym;           /* symbols of alph: "ACGT-RYMKSWHBVDN*~" for aminos */
   char     inmap[1<<8];   /* map: index -> char value */
   char     outmap[1<<8];  /* map: char value -> index */
} ALPHABET;

/* position specific state probabilities */
typedef struct {
   /* match emission probabilities for each amino acid */
   float    match[NUM_AMINO];
   /* insert emission probabilities for each amino acid */  
   float    insert[NUM_AMINO];
   /* transition state probabilities (default same as COMPO) */
   float    trans[NUM_TRANS_STATES];   
} HMM_NODE;

/* HMM Background Composition */
typedef struct {
   /* hard-coded background residue frequencies for each amino acid */
   float    freq[NUM_AMINO];
   /* background residue frequencies of the given hmm model (mean composition) */
   float    compo[NUM_AMINO];
   /* insert emission probabilities for each amino acid (uniform across positions) */
   float    insert[NUM_AMINO];
   /* transition state probabilities (default same as COMPO) */
   float    trans[NUM_TRANS_STATES];   
   /* move, loop for special transition states */
   float    spec[NUM_SPECIAL_STATES][NUM_SPECIAL_TRANS];
} HMM_COMPO;

/* HMM File Data */
typedef struct {
   int         N;       /* number of states in model */
   float*      pi;      /* initial (begin) distribution ( 0..M ) */
   float**     t;       /* state transition probabilities ( N x (N+1) ) */
   float**     e;       /* emmission probabilities ( M x K ) */

   float**     eo;      /* emission odds ratio ( M x K' ) */

   ALPHABET*   abc;     /* alphabet */           
   int         K;       /* size of alphabet */   
} HMM;

/* HMM Profile */
typedef struct {
   /* META DATA */

   /* file data */
   char*          filepath;         /* path to the file containing hmm */
   long           b_offset;         /* offset within file to beginning of hmm */
   long           e_offset;         /* offset within file to ending of hmm */

   int            numberFormat;     /* whether numbers are in real, logspace, or log-odds */

   /* entry data */
   char*          format;           /* hmm file format */
   char*          name;             /* unique name field of model in file */
   char*          acc;              /* unique accession field of model in file  */
   char*          desc;             /* description field of model in file */
   char*          alph;             /* alphabet: only amino accepted */;  

   int            nseq;             /* number of sequences used to created hmm */
   float          effn;             /* */
   long           checksum;         /* file checksum */

   float          ga[2];            /* */
   float          tc[2];            /* */
   float          nc[2];            /* */

   /* profile settings */
   int            mode;             /* enumerated search mode */
   bool           isLocal;          /* local or global? */
   bool           isMultihit;       /* multi hit or single hit? */   
   /* jump value for configuring HMM */
   int            num_J; /* number of jumps allowed by model (single hit = 1) */

   /* distribution parameters for scoring */
   DIST_PARAM     msv_dist;         /* Parameters for the Distribution for Ungapped Viterbi Scores */
   DIST_PARAM     viterbi_dist;     /* Parameters for the Distribution for Viterbi Scores */
   DIST_PARAM     forward_dist;     /* Parameters for the Distribution for Fwd/Bck Scores */

   /* MODEL DATA */
   int            N;                /* profile length (currently in use) */
   int            Nalloc;           /* profile length (allocated memory) */
   int            alph_type;        /* enumerated alphabet type */
   int            alph_leng;        /* alphabet length: AMINO = 20, DNA = 4 */
   char*          consensus;        /* consensus sequence */

   HMM_COMPO*     bg_model;         /* background composition */
   HMM_NODE*      hmm_model;        /* array of position specific probabilities */
} HMM_PROFILE;

/* Vector  */
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

/* sequence  */
typedef struct {
   int      N;          /* length of sequence */
   int      Nalloc;     /* allocated memory length */
   /* meta data */
   char*    filename;   /* filename of sequence */
   char*    name;       /* */
   char*    alph;       /* */
   char*    seq;        /* */
} SEQUENCE;

/* */
typedef struct {
   int      R;       /* number of columns = length of query */ 
   int      C;       /* number of rows = number of special states */
   int      Nalloc;  /* flat length of matrix = rows x cols */
   float*   data;    /* */
} MATRIX_2D;

/* */
typedef struct {
   int      R;       /* number of rows = length of query */
   int      C;       /* number of columns = length of target  */
   int      N;       /* number of 3rd dim = number of normal states */
   int      Nalloc;  /* number of total cells alloc'd */
   float*   data;    /* matrix cells */
} MATRIX_3D;

/* */
typedef struct {
   /* file paths */
   char*    t_filepath;          /* filepath to target (hmm or fasta) file */
   char*    q_filepath;          /* filepath to query (fasta) file */
   /* index paths */
   char*    t_indexpath;            /* index filepath for quick access of target (hmm) file */
   char*    q_indexpath;            /* index filepath for quick access of query (fasta) file */
   /* mmseqs results path */
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
   int      is_testing;             /* determines whether debug statements appear */

   /* if viterbi is precomputed, gives starting and ending coords (single result) */
   COORDS   beg;                    /* beginning coordinates of viterbi alignment */
   COORDS   end;                    /* ending coordinates of viterbi alignment */

   /* target/query metadata */
   int      t_filetype;          /* enumerated FILETYPE of target file */
   int      q_filetype;          /* enumerated FILETYPE of query file */
   /* file id number */
   int      t_fileno;            /*  */
   int      q_fileno;            /*  */ 
   /* offset into file */
   int      t_offset;
   int      q_offset;

   /* threshold scores for pipeline */
   float    viterbi_sc_threshold;
   float    fwdbck_sc_threshold;
   float    cloud_sc_threshold;
} ARGS;

/* */
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

/* */
typedef struct {
   int         id;            /* id number, determined by order in file */
   char*       name;          /* Name of HMM/FASTA in file */
   long        offset;        /* Positional offset of HMM/FASTA into file */
} F_INDEX_NODE;

/* */
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

/* */
typedef struct {
   int      num_args;      /* number of arguments */
   int      data_type;     /* data type of arguments */

   char*    name;          /* name of flag o */
   char*    long_flag;     /* long "--" flag */
   char*    short_flag;    /* single character "-" flag */
   char*    desc;          /* description of flag */
   // void*    arg_loc;       /* pointer to the location in ARGS to store option argument */
} FLAG_CMD;

/* */
typedef struct {
   int      target_id;
   int      query_id;

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

   float    cloud_fwd_sc;
} RESULT;

/* */
typedef struct {
   int      N;
   int      Nalloc;
   RESULT*  data;
   char*    filepath;
} RESULTS;

/* substitution/scoring matrix */
typedef struct {
   char*    filename;
   char*    alph;
   int      alph_len;
   int      map[256];
   float*   scores;
} SCORE_MATRIX;

/* bools of all tasks to be executed by WORKER */
/* flags set by pipeline and then passed to a generic workflow */
typedef struct {
   bool     time;             /* are we timing given tasks? */

   bool     naive_cloud;      /* */

   bool     quadratic;        /* are we running any quadratic-space algorithms? */
   bool     quad_fwdbck;      /* forward-backward */
   bool     quad_vit;         /* viterbi */
   bool     quad_trace;       /* traceback of viterbi */
   bool     quad_cloud;       /* cloud pruning search */

   bool     linear;           /* are we running any linear algorithms? */
   bool     linear_fwdbck;    /* forward-backward */
   bool     linear_vit;       /* viterbi */
   bool     linear_trace;     /* traceback of viterbi */
   bool     linear_cloud;     /* cloud pruning search */
} TASKS;

/* worker contains the necessary data structures to conduct search */
/* unnecessary components will be left NULL */
typedef struct {
   /* meta data */
   ARGS*          args;
   TASKS*         tasks;

   /* indexes of query and target data files */
   F_INDEX*       q_index;
   F_INDEX*       t_index;
   /* file pointers to query and target data file */ 
   FILE*          q_file;
   FILE*          t_file;

   /* current query and target data */
   SEQUENCE*      q_seq;
   SEQUENCE*      t_seq;
   HMM_PROFILE*   t_prof;

   /* edgebounds for cloud search (linear space) */
   EDGEBOUNDS*    edg_fwd;
   EDGEBOUNDS*    edg_bck;
   EDGEBOUNDS*    edg_diag;
   EDGEBOUNDS*    edg_row;

   /* alignment traceback for viterbi */
   ALIGNMENT*     traceback;

   /* dynamic programming matrices */
   MATRIX_3D*     st_MX;  /* normal state matrix (quadratic space) */
   MATRIX_3D*     st_MX3; /* normal state matrix (linear space) */
   MATRIX_2D*     sp_MX;  /* special state matrix (quadratic space) */

   /* times for tasks */
   TIMES*         times;
   /* scores for algorithms */
   SCORES*        scores;
   /* results for worker */
   RESULTS*       results;
   /* clock for taking times */
   CLOCK*         clock;
} WORKER;


/* === GLOBAL VARIABLES === */

extern char*   STATE_NAMES[];
extern char*   STATE_FULL_NAMES[];

extern char*   PIPELINE_NAMES[];
extern void    (*PIPELINES[])(WORKER*);

extern char*   MODE_NAMES[];
extern char*   VERBOSITY_NAMES[];

extern char*   ALPHABET_NAMES[];
extern int     ALPHABET_LENGTHS[];

extern char*   FILE_TYPE_EXTS[];
extern char*   FILE_TYPE_NAMES[];
extern int     FILE_TYPE_MAP[];

/* alphabetically-ordered amino lookup and reverse lookup tables */
extern char    ALPH_AMINO_CHARS[];
extern char    AA[];
extern int     AA_REV[];
/* background frequencies of null model, normal and log space */
extern double  BG_MODEL[];
extern double  BG_MODEL_log[];

#endif /* _STRUCTS_H */
