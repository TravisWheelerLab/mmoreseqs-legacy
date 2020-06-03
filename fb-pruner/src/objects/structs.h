/*******************************************************************************
 *     FILE:   structs.h
 *  PURPOSE:   All Data Structures used by Cloud Search.
 *
 *  AUTHOR:    Dave Rich
 *     BUG:    
 *******************************************************************************/

#ifndef _STRUCTS_H
#define _STRUCTS_H

/* === STDLIB DATA TYPES === */
#include <stdbool.h>
#include <time.h>
#include <sys/types.h>

/* === MACRO FUNCTIONS === */
#include "objects/structs_macros.h"

/* === ENUMERATIONS === */
#include "objects/structs_enums.h"

/* === DECLARE HIDDEN FUNCTIONS === */
/* (*** c99 standard hides this function from <string.h> ***) */
// extern char* strdup(const char*);
// extern ssize_t getline(char **lineptr, size_t *n, FILE *stream);

/* === STRUCTS === */

/* integer ranges */
typedef struct {
   int beg;
   int end;
} RANGE;

/* coordinates in matrix */
typedef struct {
   int i; /* row index */
   int j; /* col index */
} COORDS;

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
   int            edg_mode;   /* whether edges are stored by-row or by-diag */

   /* dimensions of embedding matrix */
   int            Q;          
   int            T;
} EDGEBOUNDS;

/* a vector of each set of bounds for cloud search space which account a specific row */
typedef struct {
   int            N;          /* current size of array */
   int            Nalloc;     /* allocated size of array */
   int*           rows_N;     /* current number of bounds in row */
   BOUND*         rows;       /* array of bounds, each for a specific row */
   int            row_max;    /* maximum number of bounds in row */

   /* dimension of embedding matrix */
   int            Q;
   int            T;
} EDGEBOUND_ROWS;

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

   /* dimensions of embedded matrix */
   int         Q;
   int         T;
} ALIGNMENT;

/* clock for timing events (wrapper for squid stopwatch above) */
typedef struct {
   double   start;            /* captures start time */
   double   stop;             /* captures stop time */
   double   duration;         /* captures difference between start and stop */

   int      N;                /* current utilized stamp length*/
   int      Nalloc;           /* allocated stamp length */
   double*  stamps;           /* for storing multiple durations */
} CLOCK;

/* distribution parameters */
typedef struct {
   float    param1;
   float    param2;
} DIST_PARAM;

/* alphabet (amino, dna, etc) (modeled after EASEL) */
typedef struct {
   int      K;             /* size of unique alphabet */
   int      Kp;            /* size of total symbols: alphabet + special symbols  */
   char*    sym;           /* symbols of alph: "ACGT-RYMKSWHBVDN*~" for aminos */
   char     inmap[1 << 8]; /* map: index -> char value */
   char     outmap[1 << 8]; /* map: char value -> index */
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

   /* hmm model/profile */
   int            N;                /* profile length (currently in use) */
   int            Nalloc;           /* profile length (allocated memory) */
   int            alph_type;        /* enumerated alphabet type */
   int            alph_leng;        /* alphabet length: AMINO = 20, DNA = 4 */
   char*          consensus;        /* consensus sequence */

   HMM_COMPO*     bg_model;         /* background composition */
   HMM_NODE*      hmm_model;        /* array of position specific probabilities */
} HMM_PROFILE;

/* vector of edgebounds  */
typedef struct {
   BOUND*   data;    /* array of data type */
   int      N;       /* current length of array in use */
   int      Nalloc;  /* current length of array allocated */
}  VECTOR_BOUND;

/* vector of structs */
typedef struct {
   TRACE*   data;     /* array of data type */
   int      N;        /* current length of array in use */
   int      Nalloc;   /* current length of array allocated */
}  VECTOR_TRACE;

/* sequence */
typedef struct {
   int      N;          /* length of sequence */
   int      Nalloc;     /* allocated memory length */
   /* meta data */
   char*    filename;   /* filename of sequence */
   char*    name;       /* */
   char*    alph;       /* */
   char*    seq;        /* */
} SEQUENCE;

/* 2-dimensional float matrix */
typedef struct {
   int      R;       /* number of columns = length of query */
   int      C;       /* number of rows = number of special states */
   int      Nalloc;  /* flat length of matrix = rows x cols */
   float*   data;    /* */
   bool     clean;   /* whether data has been cleared / all cells set to -INF */
} MATRIX_2D;

/* 3-dimensional float matrix */
typedef struct {
   int      R;       /* number of rows = length of query */
   int      C;       /* number of columns = length of target  */
   int      N;       /* number of 3rd dim = number of normal states */
   int      Nalloc;  /* number of total cells alloc'd */
   float*   data;    /* matrix cells */
   bool     clean;   /* whether data has been cleared / all cells set to -INF */
} MATRIX_3D;

/* commandline arguments */
typedef struct {
   bool     is_index_given;         /* is an index file supplied? */

   /* file paths */
   char*    t_filepath;             /* filepath to target (hmm or fasta) file */
   char*    q_filepath;             /* filepath to query (fasta) file */
   /* index paths */
   char*    t_indexpath;            /* index filepath for quick access of target (hmm) file */
   char*    q_indexpath;            /* index filepath for quick access of query (fasta) file */
   /* output path */
   char*    output_filepath;        /* filepath to output results to; "!stdout" => stdout */
   /* target/query metadata */
   int      t_filetype;             /* enumerated FILETYPE of target file */
   int      q_filetype;             /* enumerated FILETYPE of query file */
   /* for specified range of targets/queries in file */
   RANGE    t_range;                /* start-end range of targets in file (inclusive) */
   RANGE    q_range;                /* start-end range of queries in file (inclusive) */
   /* temporary work folder */
   char*    tmp_folderpath;         /* location to build a temporary work folder */

   /* mmseqs-plus search params */
   char*    mmseqs_res_filepath;       /* filepath to mmseqs .m8 results file */
   char*    mmseqs_plus_filepath;      /* filepath to mmseqs .m8+ results file */
   char*    mmseqs_tmp_filepath;       /* filepath to mmseqs temporary files root */
   /* mmseqs reference files */
   char*    t_lookup_filepath;      /* filepath to mmseqs target lookup file */
   char*    q_lookup_filepath;      /* filepath to mmseqs query lookup file */
   /* mmseqs specified result file range */
   RANGE    mmseqs_range;

   /* cloud search tuning vars */
   float    alpha;                  /* cloud search: x-drop pruning ratio */
   float    alpha_max;              /* cloud search: x-drop maximum drop before termination */
   int      beta;                   /* cloud search: number of antidiag passes before pruning  */

   /* pipeline options */
   int      pipeline_mode;          /* which workflow pipeline to use */
   int      verbose_level;         /* levels of verbosity */
   int      search_mode;            /* alignment search mode */
   bool     is_testing;             /* determines whether debug statements appear */

   /* if viterbi is precomputed, gives starting and ending coords (single result) */
   COORDS   beg;                    /* beginning coordinates of viterbi alignment */
   COORDS   end;                    /* ending coordinates of viterbi alignment */

   /* threshold scores for pipeline */
   float    viterbi_threshold;
   float    fwdbck_threshold;
   float    cloud_threshold;
} ARGS;

/* scores */
typedef struct {
   /* naive algs */
   float    naive_cloud_fwd;
   float    naive_cloud_bck;
   /* quadratic algs */
   float    quad_vit;
   float    quad_trace;
   float    quad_fwd;
   float    quad_bck;
   float    quad_cloud_fwd;
   float    quad_cloud_bck;
   /* linear */
   float    lin_vit;
   float    lin_trace;
   float    lin_fwd;
   float    lin_bck;
   float    lin_cloud_fwd;
   float    lin_cloud_bck;
   /* statistics */
   float    perc_cells;
   float    perc_window;
} SCORES;

/* times to execute given operations */
typedef struct {
   /* load times */
   double    load_target_index;
   double    load_query_index;

   double    load_target;
   double    load_query;

   /* linear algs */
   double    lin_vit;
   double    lin_trace;
   double    lin_fwd;
   double    lin_bck;
   double    lin_cloud_fwd;
   double    lin_cloud_bck;
   double    lin_merge;
   double    lin_reorient;
   double    lin_bound_fwd;
   double    lin_bound_bck;
   double    lin_total_cloud;

   /* quadratic algs */
   double    quad_vit;
   double    quad_trace;
   double    quad_fwd;
   double    quad_bck;
   double    quad_cloud_fwd;
   double    quad_cloud_bck;
   double    quad_merge;
   double    quad_reorient;
   double    quad_bound_fwd;
   double    quad_bound_bck;  
   double    quad_total_cloud;

   /* naive algs */
   double    naive_bound_fwd;
   double    naive_bound_bck;  
} TIMES;

/* times to execute given operations */
typedef struct {
   /* load times */
   long    load_target_index;
   long    load_query_index;

   long    load_target;
   long    load_query;

   /* linear algs */
   long    lin_vit;
   long    lin_trace;
   long    lin_fwd;
   long    lin_bck;
   long    lin_cloud_fwd;
   long    lin_cloud_bck;
   long    lin_merge;
   long    lin_reorient;
   long    lin_bound_fwd;
   long    lin_bound_bck;

   /* quadratic algs */
   long    quad_vit;
   long    quad_trace;
   long    quad_fwd;
   long    quad_bck;
   long    quad_cloud_fwd;
   long    quad_cloud_bck;
   long    quad_merge;
   long    quad_reorient;
   long    quad_bound_fwd;
   long    quad_bound_bck;  

   /* naive algs */
   long    naive_bound_fwd;
   long    naive_bound_bck;  
} TIMES_RAW;

/* hmm location within file */
typedef struct {
   int         id;            /* id number, determined by order in file */
   char*       name;          /* Name of HMM/FASTA in file */
   long        offset;        /* Positional offset of HMM/FASTA into file */
   int         mmseqs_id;     /* id number, referencing mmseqs lookup file */
} F_INDEX_NODE;

/* index for offset locations into a file, searchable by name or id */
typedef struct {
   int            N;          /* Number of location index nodes used */
   int            Nalloc;     /* Number of location index nodes allocated */
   F_INDEX_NODE*  nodes;      /* List of nodes for location of each HMM/FASTA in file */    

   char*          index_path;    /* Filepath of index file (NULL if built and not loaded) */
   char*          lookup_path;   /* Filepath to mmseqs lookup file (NULL if not used) */
   char*          source_path;   /* Filepath of file being indexed */
   char*          delim;         /* one of more delimiter of header fields */

   int            filetype;      /* Type of file being indexed (HMM, FASTA, etc) */
   int            sort_type;     /* Whether the index nodes list has been sorted, and by which field */
   int            mmseqs_names;  /* Whether index is using names from mmseqs lookup */
} F_INDEX;


/* flags for command line arguments */
typedef struct {
   char*    name;          /* name of flag */

   int      num_args;      /* number of arguments */
   int      data_type;     /* data type of arguments */
   void*    arg_loc;       /* pointer to the location in ARGS to store option argument */

   char*    long_flag;     /* long "--" flag */
   char*    short_flag;    /* single character "-" flag */
   char*    desc;          /* description of flag */
} FLAG_CMD;

/* results fields */
typedef struct {
   /* result unique identifier */
   int      result_id;
   /* cloud search lookup index id */
   int      target_id;
   int      query_id;
   /* mmseqs lookup index id */
   int      target_mid;
   int      query_mid;
   /* name */
   char*    target_name;
   char*    query_name;
   /* */
   float    perc_id;
   int      aln_len;
   int      mismatch;
   int      gap_openings;
   /* alignment window */
   COORDS   beg;
   COORDS   end;
   /* alignment window */
   int      query_start;
   int      query_end;
   int      target_start;
   int      target_end;
   /* */
   double   e_value;
   int      bit_score;
   /* */
   float    cloud_fwd_sc;
   /* number of cells computed */
   int      cpu_cloud_cells;
   int      cloud_cells;
   int      cpu_total_cells;
   int      total_cells;
} RESULT;

/* array of results */
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

/* data struct for passing debugger data */
typedef struct {
   /* debug options */
   bool           is_debugging;     /* debug flag determines whether testing debug output is generated */
   int            verbose_level;   /* level of verbosity; amount to output */  
   bool           is_embed;         /* whether to embed linear space matrices in quadratic matrices */
   bool           is_viz;           /* whether to fill data viz matrix */
   /* output paths */
   char*          dbfp_path;        /* debugger ouput filepath */
   FILE*          dbfp;             /* debugger output file pointer */
   /* data structs */
   MATRIX_3D*     test_MX;          /* quadratic testing matrix */
   MATRIX_2D*     cloud_MX;         /* cloud matrix for visualizing */
} DEBUG_KIT;

/* bools of all tasks to be executed by WORKER */
/* flags set by pipeline and then passed to a generic workflow */
typedef struct {
   /* results output */
   bool     scores;           /* are we reporting scores for given tasks? */
   bool     time;             /* are we reporting times for given tasks? */

   /* output heatmaps */
   bool     heatmaps;         /* output heatmaps */

   /* mmseqs lookup */
   bool     mmseqs_lookup;    /* if we have a mmseqs name lookup table */

   /* naive algs */
   bool     naive_cloud;      /* naive cloud search */

   /* quadratic algs */
   bool     quadratic;        /* are we running any quadratic-space algorithms? */
   bool     quad_fwd;         /* forward */
   bool     quad_bck;         /* backward */
   bool     quad_vit;         /* viterbi */
   bool     quad_trace;       /* traceback of viterbi */
   bool     quad_bound_fwd;   /* bounded forward (requires cloud) */
   bool     quad_bound_bck;   /* bounded backward (requires cloud) */

   /* linear algs */
   bool     linear;           /* are we running any linear algorithms? */
   bool     lin_fwd;          /* forward-backward */
   bool     lin_bck;          /* backward */
   bool     lin_vit;          /* viterbi */
   bool     lin_trace;        /* traceback of viterbi */
   bool     lin_bound_fwd;    /* bounded forwarded (requires cloud) */
   bool     lin_bound_bck;    /* bounded backward (requires cloud) */
} TASKS;

/* bools of which scores and times to be reported */
typedef struct {
   /* TIMES */
   /* linear algs */
   bool     lin_fwd_t;
   bool     lin_bck_t;
   bool     lin_vit_t;
   bool     lin_trace_t;
   bool     lin_cloud_fwd_t;
   bool     lin_cloud_bck_t;
   bool     lin_merge_t;
   bool     lin_reorient_t;
   bool     lin_bound_fwd_t;
   bool     lin_bound_bck_t;
   /* quadratic algs */
   bool     quad_fwd_t;
   bool     quad_bck_t;
   bool     quad_vit_t;
   bool     quad_trace_t;
   bool     quad_cloud_fwd_t;
   bool     quad_cloud_bck_t;
   bool     quad_merge_t;
   bool     quad_reorient_t;
   bool     quad_bound_fwd_t;
   bool     quad_bound_bck_t;

   /* SCORES */
   /* linear algs */
   bool     lin_fwd_sc;
   bool     lin_bck_sc;
   bool     lin_vit_sc;
   bool     lin_bound_fwd_sc;
   bool     lin_bound_bck_sc;  
   /* quadratic algs */
   bool     quad_fwd_sc;
   bool     quad_bck_sc;
   bool     quad_vit_sc;
   bool     quad_bound_fwd_sc;
   bool     quad_bound_bck_sc;  
} REPORT;

/* collection of all tuning parameters for cloud search */
typedef struct {
   float       alpha;         /* x-drop for local,  */
   float       alpha_max;     /* x-drop for global, determines when to terminate search. looser (larger) than alpha. */
   float       beta;          /* number of traversed antidiags before pruning begins */
} CLOUD_PARAMS;

/* worker contains the necessary data structures to conduct search */
/* unnecessary components will be left NULL */
typedef struct {
   /* meta data */
   ARGS*          args;
   TASKS*         tasks;
   /* fields to report */
   REPORT*        report;
   /* output file pointer */
   FILE*          out_file;

   /* indexes of query and target data files */
   F_INDEX*       q_index;
   F_INDEX*       t_index;
   /* file pointers to query and target data file */
   FILE*          q_file;
   FILE*          t_file;
   /* id of currently loaded query and target */
   int            q_id;
   int            t_id;

   /* current query and target data */
   SEQUENCE*      q_seq;
   SEQUENCE*      t_seq;
   HMM_PROFILE*   t_prof;

   /* edgebounds for cloud search */
   EDGEBOUNDS*    edg_fwd;
   EDGEBOUNDS*    edg_bck;
   EDGEBOUNDS*    edg_diag;
   EDGEBOUNDS*    edg_row;

   /* edgebound row object; helper for reorientation of edgebounds */
   EDGEBOUND_ROWS*   edg_rows_tmp;

   /* cloud pruning parameters */
   CLOUD_PARAMS  cloud_params;

   /* alignment traceback for viterbi */
   ALIGNMENT*     traceback;

   /* dynamic programming matrices */
   MATRIX_3D*     st_MX;         /* normal state matrix (quadratic space) */
   MATRIX_3D*     st_MX3;        /* normal state matrix (linear space) */
   MATRIX_2D*     sp_MX;         /* special state matrix (quadratic space) */
   MATRIX_3D*     st_cloud_MX;   /* matrix for naive cloud search */

   /* times for tasks */
   TIMES*         times;
   TIMES_RAW*     times_raw;
   /* scores for algorithms */
   SCORES*        scores;
   /* results from mmseqs */
   RESULTS*       results_in;
   /* results to output */
   RESULTS*       results; 
   /* current result */
   RESULT*        result;
   /* clok for taking times */
   CLOCK*         clok;
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

/* commandline arg objects */
extern ARGS*      args;
extern int        num_flag_cmds;
extern FLAG_CMD   COMMAND_OPTS[];
extern char*      DATATYPE_NAMES[];

/* debugging data */
extern DEBUG_KIT* debugger;

#endif /* _STRUCTS_H */
