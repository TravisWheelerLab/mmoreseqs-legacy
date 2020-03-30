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
#include "objects/structs_macros.h"

/* === ENUMERATIONS === */
#include "objects/structs_enums.h"

/* === STRUCTS === */

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

typedef struct {
   float    param1;
   float    param2;
} DIST_PARAM;

/* position specific state probabilities */
typedef struct {
   /* match emission probabilities for each amino acid */
   float    match[NUM_AMINO];
   /* insert emission probabilities for each amino acid */  
   float    insert[NUM_AMINO];
   /* transition state probabilities (default same as COMPO) */
   float    trans[NUM_TRANS_STATES];   
} HMM_NODE;

/* HMM Background Probabilities */
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
   int      num_J; /* number of jumps allowed by model */
} HMM_BG;

/* HMM File Data */
typedef struct {
   
} HMM;

/* HMM Profile */
typedef struct {
   /* META DATA */
   /* file data */
   char*          filepath;         /* path to the file containing hmm */
   int            b_offset;         /* offset within file to beginning of hmm */
   int            e_offset;         /* offset within file to ending of hmm */

   /* entry data */
   char*          name;             /* unique name field of model in file */
   char*          acc;              /* unique accession field of model in file  */
   char*          desc;             /* description field of model in file */
   char*          alph;             /* PROTEIN or DNA (Only PROTEIN ACCEPTED ) */;  

   /* profile settings */
   int            mode;             /* enumerated mode */
   int            isLocal;          /* local or global? */
   int            isMultihit;       /* multi hit or single hit? */   

   /* distribution parameters for scoring */
   DIST_PARAM*    msv_dist;         /* Parameters for the Distribution for Ungapped Viterbi Scores */
   DIST_PARAM*    viterbi_dist;     /* Parameters for the Distribution for Viterbi Scores */
   DIST_PARAM*    forward_dist;     /* Parameters for the Distribution for Fwd/Bck Scores */

   /* MODEL DATA */
   int            N;                /* profile length (currently in use) */
   int            Nalloc;           /* profile length (allocated memory) */
   int            alph_leng;        /* alphabet length: AMINO = 20, DNA = 4 */

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

   /* threshold scores for pipeline */
   float    viterbi_sc_threshold;
   float    fwdbck_sc_threshold;
   float    cloud_sc_threshold;
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
   int         id;            /* id number, determined by order in file */
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
   int      data_type;     /* data type of arguments */

   char*    name;          /* name of flag o */
   char*    long_flag;     /* long "--" flag */
   char*    short_flag;    /* single character "-" flag */
   char*    desc;          /* description of flag */
   // void*    arg_loc;       /* pointer to the location in ARGS to store option argument */
} FLAG_CMD;


typedef struct {
   int      target_int;
   int      query_int;

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

extern char*   FILE_TYPE_EXTS[];
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
