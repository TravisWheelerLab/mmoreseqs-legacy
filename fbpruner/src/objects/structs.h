/*******************************************************************************
 *     FILE:   structs.h
 *  PURPOSE:   All Data Structures used by Cloud Search.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _STRUCTS_H
#define _STRUCTS_H

/* === STDLIB DATA TYPES === */
#include <stdbool.h>
#include <time.h>
#include <sys/types.h>

/* === MACRO FUNCTIONS === */
#include "structs_macros.h"

/* === ENUMERATIONS === */
#include "structs_enums.h"

/* === DECLARE HIDDEN FUNCTIONS (for older compilers) === */
/* (*** c99 standard hides this function from <string.h> ***) */
// extern char* strdup(const char*);
// extern ssize_t getline(char **lineptr, size_t *n, FILE *stream);

/* === EASEL === */
#include "easel.h"
#include "esl_hmm.h"
#include "esl_alphabet.h"
#include "esl_sq.h"

/* === STRUCTS === */

/* template datatype */
#define XXX          int
/* alias primitive datatypes */
#define INT          int
#define FLT          float
#define CHAR         char 
#define DBL          double
#define STR          char*
/* datatype for matrices */
#define DATA         float

/* generic datatype that can hold most basic datatypes */
typedef union {
   int      i;
   float    f;
   double   d;
   long     l;
   char*    s;
   char     c;
   bool     b;
} GENERIC;

/* generic data with specified type */
typedef struct {
   GENERIC        data;
   DATATYPES      type;
} GEN_DATA;

/* coordinates in matrix */
typedef struct {
   int      i;       /* row index */
   int      j;       /* col index */
} COORDS;

/* integer ranges */
typedef struct {
   int      beg;      /* begin: left bound */
   int      end;      /* end: right bound */
} RANGE;

/* bound for single row or diagonal */
typedef struct {
   int      id;         /* current anti-diagonal OR row */
   int      lb;         /* bottom-left (lower) edge-bound range */
   int      rb;         /* top-right (upper) edge-bound range */
} BOUND;

/* string object */
typedef struct {
   int      N;
   char*    data;
} X_STRING;

/* given cell of alignment */
/* TODO: change to t_0, q_0 indexing */
typedef struct {
   int      q_0;       /* index in query */
   int      t_0;       /* index in target */
   int      st;      /* state at index */
} TRACE;

/* === VECTORS === */

/* default size of vectors when created */
#define VECTOR_INIT_SIZE 32

/* vector of template data */
typedef struct {
   XXX*        data;     /* array of data type */
   int         N;        /* current length of array in use */
   int         Nalloc;   /* current length of array allocated */
}  VECTOR_XXX;

/* vector of characters */
typedef struct {
   CHAR*    data;             /* array of data type */
   int      N;                /* current length of array in use */
   int      Nalloc;           /* current length of array allocated */
}  VECTOR_CHAR;

/* vector of strings */
typedef struct {
   STR*     data;             /* array of data type */
   int      N;                /* current length of array in use */
   int      Nalloc;           /* current length of array allocated */
}  VECTOR_STR;


/* vector of integers */
typedef struct {
   INT*     data;             /* array of data type */
   int      N;                /* current length of array in use */
   int      Nalloc;           /* current length of array allocated */
}  VECTOR_INT;

/* vector of floats */
typedef struct {
   FLT*     data;             /* array of data type */
   int      N;                /* current length of array in use */
   int      Nalloc;           /* current length of array allocated */
}  VECTOR_FLT;

/* vector of doubles */
typedef struct {
   DBL*     data;             /* array of data type */
   int      N;                /* current length of array in use */
   int      Nalloc;           /* current length of array allocated */
}  VECTOR_DBL;

/* vector of edgebounds  */
typedef struct {
   BOUND*   data;    /* array of data type */
   int      N;       /* current length of array in use */
   int      Nalloc;  /* current length of array allocated */
}  VECTOR_BOUND;

/* vector of range structs */
typedef struct {
   RANGE*   data;     /* array of data type */
   int      N;        /* current length of array in use */
   int      Nalloc;   /* current length of array allocated */
}  VECTOR_RANGE;

/* vector of trace structs */
typedef struct {
   TRACE*   data;     /* array of data type */
   int      N;        /* current length of array in use */
   int      Nalloc;   /* current length of array allocated */
}  VECTOR_TRACE;

/* vector of data structs */
typedef struct {
   DATA*    data;     /* array of data type */
   int      N;        /* current length of array in use */
   int      Nalloc;   /* current length of array allocated */
}  VECTOR_DATA;

/* === MAP === */

/* node for UMAP */
typedef struct {
   char*    key;
   XXX      value;
   int      type;
} XXX_UMAP_NODE;

/* unordered map */
typedef struct {
   int               N;       /* number of nodes */
   int               Nalloc;  /* number of nodes allocated */
   XXX_UMAP_NODE*    nodes;   /* key-value pairs */
} XXX_UMAP;

/* reader for parsing files */
typedef struct {
   FILE*             fp;            /* file pointer to be read */
   char*             filename;      /* string to file location */
   /* buffer data */
   char*             buffer;        /* character buffer */
   size_t            N;             /* length of occupied buuffer */
   size_t            Nalloc;        /* length of allocated buffer */
   int               lines_read;    /* number of lines read */
   char*             token;         /* pointer to current token in tokenized buffer */
   VECTOR_CHAR*      tokens;        /* pointers to all tokens in buffer line */
   char*             delim;         /* current delimiter being used for parsing tokens */
   /* location in file */
   long int          cur_offset;    /* current offset of start of buffer into file */
   long int          file_size;     /* size of file (in chars) */
   bool              is_eof;        /* is file at end of file */
} READER;

/* single option */
typedef struct {
   char*             opt_long;   /* long option */
   char*             opt_short;  /* short option */
   int               arg_type;   /* datatype of arguments */
   int               arg_num;    /* number of args allowed */
   char*             help_info;  /* report help data */
} OPT;

/* all command line options */
typedef struct {
   /* option vector */
   OPT*              data;
   int               N;
   int               Nalloc;
   char*             help_info;
} OPTS;

/* set of bounds for cloud search space */
typedef struct {
   /* size of edgebounds */
   int            N;             /* current size of bounds array */
   int            Nalloc;        /* allocated size of bounds array */
   /* dimensions of embedding matrix */
   int            Q;          
   int            T;
   /* whether edges are stored by-row or by-diag */
   int            edg_mode;   
   /* row/antidiagonal indexes */
   VECTOR_INT*    ids;           /* identity of each row/diag */
   VECTOR_INT*    ids_idx;       /* array of indexes the heads of each row/diag */
   /* data */
   BOUND*         bounds;        /* array of bounded ranges along a row/diag */
} EDGEBOUNDS;

/* a vector of each set of bounds for cloud search space which account a specific row */
typedef struct {
   /* size */
   int            N;          /* current size of array */
   int            Nalloc;     /* allocated size of array */
   int*           rows_N;     /* current number of bounds in row */
   int            row_max;    /* maximum number of bounds in row */
   /* dimension of embedding matrix */
   int            Q;
   int            T;
   /* data */
   BOUND*         rows;       /* array of bounds, each for a specific row */
} EDGEBOUND_ROWS;

/* alignment for viterbi traceback */
typedef struct {
   /* consensus data (for alignment output) */
   VECTOR_INT*    seq_beg;       /* index of beginnings of each alignment in sequence */
   VECTOR_INT*    seq_end;       /* index of endings of each alignment in sequence */
   VECTOR_CHAR*   sequence;      /* character sequence showing alignment between two sequences */
   /* data */
   VECTOR_INT*    tr_beg;        /* index of every alignment begin (B state) in traceback */
   VECTOR_INT*    tr_end;        /* index of every alignment end (E state) in traceback */
   VECTOR_TRACE*  traces;        /* list of all (state,i,j) TRACES in ALIGNMENT */
   /* mmseqs' cigar-style alignments */
   char*          cigar_aln;     /* ex: 379M173M2D41M2D65M6I21 (see MMseqs userguide pg.35) */
   /* hmmer-style alignments */
   char*          target_aln;    /* |qdrfLePqcrilDlevWdqCYfRWlPvLeikgGGqpq| */
   char*          center_aln;    /* |++  L    ++  l+vW qCYfRW+P  e+ +G    | */
   char*          query_aln;     /* |DEDLLSMDATVKALSVWTQCYFRWVPKAEVVNGDPGT| */
   char*          state_aln;     /* |MMMMIIIDIDMM...                      | */
   /* trace endpoints */
   int            beg;           /* current beginning index in traces */
   int            end;           /* current end index in traces */
   /* dimensions of embedded matrix */
   int            Q;
   int            T;
   /* counts of hits, mismatches, and gaps */
   int            num_gaps;
   int            num_misses;
   int            num_matches;
   float          perc_id;
   /* meta data */
   int            aln_len;       /* alignment length */
} ALIGNMENT;

/* clock for timing events (wrapper for squid stopwatch above) */
typedef struct {
   /* current time */
   double      start;            /* captures start time */
   double      stop;             /* captures stop time */
   double      duration;         /* captures difference between start and stop */
   /* previous times */
   double      program_start;    /* time captured when clock is created */
   VECTOR_DBL* stamps;
} CLOCK;

/* distribution parameters */
typedef struct {
   float    param1;
   float    param2;
} DIST_PARAM;

/* alphabet (amino, dna, etc) (modeled after EASEL) */
typedef struct {
   int      K;                   /* size of unique alphabet */
   int      Kp;                  /* size of total symbols: alphabet + special symbols  */
   char*    sym;                 /* symbols of alph: "ACGT-RYMKSWHBVDN*~" for aminos */
   char     inmap[(1 << 8)];     /* map: index -> char value */
   char     outmap[(1 << 8)];    /* map: char value -> index */
} ALPHABET;

/* Null Model for computing Composition Bias (Modeled after Easel) */
typedef struct {
   float*         f;       /* null_1 background residue frequencies [0..K-1]: set at initialization */
   float          p1;      /* null1's transition prob:  */
   ESL_HMM*       fhmm;    /* bias filter: p7_bg_SetFilter() sets this, from model's mean composition */
   float          omega;   /* the "prior" on null2/null3: set at initialization (one omega for both null types)  */
   ESL_ALPHABET*  abc;     /* reference to alphabet in use: set at initialization */
   /* NOTE: should move this to SEQUENCE? */
   ESL_SQ*        sq;      /* sequence  */
} HMM_BG;

/* position-specific state probabilities */
typedef struct {
   /* match emission probabilities for each amino acid at position */
   float    match[NUM_AMINO_PLUS_SPEC];
   /* insert emission probabilities for each amino acid at position */
   float    insert[NUM_AMINO_PLUS_SPEC];
   /* transition state probabilities at position (default same as COMPO) */
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
   /* composition for null model */
   HMM_BG*  hmm_bg;
} HMM_COMPO;

/* HMM Profile */
typedef struct {
   /* file data */
   char*          filepath;         /* path to the file containing hmm */
   long           b_offset;         /* offset within file to beginning of hmm */
   long           e_offset;         /* offset within file to ending of hmm */
   /* meta data */
   int            numberFormat;     /* whether numbers are in real, logspace, or log-odds */
   /* entry data */
   char*          format;           /* hmm file format */
   char*          name;             /* unique name field of model in file */
   char*          acc;              /* unique accession field of model in file  */
   char*          desc;             /* description field of model in file */
   char*          alph;             /* alphabet: only amino accepted */
   /* */
   int            nseq;             /* number of sequences used to created hmm */
   float          effn;             /* */
   long           checksum;         /* file checksum */
   /* cutoff thresholds */
   float          ga[2];            /* "gathering" cutoff :: threshold for family membership */
   float          nc[2];            /* "noise" cutoff :: threshold for highest known false positive */
   float          tc[2];            /* "trusted" cutoff :: threshold for lowest scoring true positive that is above all know false positives */
   /* profile settings */
   int            mode;             /* enumerated search mode */
   bool           isLocal;          /* local or global? */
   bool           isMultihit;       /* multi hit or single hit? */
   /* jump value for configuring HMM */
   float          num_J;            /* number of jumps allowed by model (single hit = 1) */
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
   /* main model */
   HMM_COMPO*     bg_model;         /* background composition */
   HMM_NODE*      hmm_model;        /* array of position specific probabilities */
   /* submodel */
   int            N_full;           /* profile length of full model */
   HMM_NODE*      hmm_model_full;   /* array of position at start of full model */
} HMM_PROFILE;

/* sequence */
typedef struct {
   int            N;          /* length of sequence (can be the length of a subsequence) */
   int            Nalloc;     /* allocated memory length */
   int            full_N;     /* length of entire sequence */
   char*          full_seq;   /* entire sequence */
   char*          seq;        /* genomic sequence (can point to the start of a subsequence) */
   int*           dseq;       /* digitized genomic sequence */
   /* imported from easel */
   ESL_SQ*        esl_dsq;    /* easel's digitized sequence */
   /* meta data */
   char*          filename;   /* filename of sequence */
   char*          header;     /* full descriptor line */
   char*          name;       /* name of sequence */
   char*          acc;        /* accession of sequence */
   char*          alph;       /* alphabet (currently only supports AMINO) */
} SEQUENCE;

/* 2-dimensional float matrix */
typedef struct {
   /* dimensions */
   int         R;       /* number of columns = length of query */
   int         C;       /* number of rows = number of special states */
   /* allocated size */
   int         Nalloc;  /* flat length of matrix = rows x cols */
   float*      data;    /* */
   bool        clean;   /* whether data has been cleared / all cells set to -INF */
} MATRIX_2D;

/* 3-dimensional float matrix */
typedef struct {
   /* dimensions */
   int         R;       /* number of 1st dim = length of query */
   int         C;       /* number of 2nd dim = length of target  */
   int         N;       /* number of 3rd dim = number of normal states */
   /* allocated size */
   int         Nalloc;  /* number of total cells alloc'd (flat dimension) */
   float*      data;    /* matrix cells */
   bool        clean;   /* whether data has been cleared / all cells set to -INF */
} MATRIX_3D;

/** 3-dimensional sparse float matrix */
typedef struct {
   /* dimensions */
   int            D1;            /* number of rows = length of query */
   int            D2;            /* number of columns = length of target  */
   int            D3;            /* number of 3rd dim = number of normal states */
   /* flat dimension and allocated size */
   int            N;             /* number of cells in sparse matrix */
   int            Nalloc;        /* number of total cells alloc'd */
   /* definitions of inner and outer edge boundaries */
   EDGEBOUNDS*    edg_inner;     /* edgebounds which describe of active cells of sparse matrix */
   EDGEBOUNDS*    edg_outer;     /* edgebounds which describe the shape of sparse matrix */
   /* offset to starts of each inner edgebound ranges (plus upper and lower rows) */
   VECTOR_INT*    imap_prv;      /* maps edg_inner to offsets into data, at previous row */
   VECTOR_INT*    imap_cur;      /* maps edg_inner to offsets into data, at current row */
   VECTOR_INT*    imap_nxt;      /* maps edg_inner to offsets into data, at next row */
   /* offset to starts of each outer edgebound ranges */
   VECTOR_INT*    omap_cur;      /* maps edg_outer to offsets into data */
   /* cell data */
   VECTOR_FLT*    data;          /* matrix cells */
   bool           clean;         /* whether data has been cleared / all cells set to -INF */
   /* iterators for traversing matrix */
   int            q_iter;        /* current query position iterator */
   int            t_iter;        /* current target position iterator */
   int            r_iter;        /* current edgebound iterator for retrieving next row range */
} MATRIX_3D_SPARSE;

/** full  */
typedef struct {
   /* dimensions */
   int                  Q;             /* query length */
   int                  T;             /* target length */
   /* normal state data {MID} */ 
   MATRIX_3D*           st_MX;         /* quadratic space normal state matrix */
   MATRIX_3D*           st_MX3;        /* linear space normal state matrix */
   MATRIX_3D_SPARSE*    st_SMX;        /* sparse space normal state matrix */
   /* special state data {ENJCB} */
   MATRIX_2D*           sp_MX;         /* special state matrix */     
   /* normal matrix type */
   bool                 is_quad;
   bool                 is_lin;
   bool                 is_sparse;
   /* meta data */
   int                  data_type;     /* whether dp is in log-scale, normal-scale, etc */
   bool                 clean;         /* if matrix is clean (all cells set to -INF)  */
   float                scaled;        /* if all cells are scaled by a common factor */                
} DP_MATRIX;

/* commandline arguments */
typedef struct {
   /* --- FULL CMD LINE --- */
   /* all commandline and opts */
   char*    cmdline;                /* full commandline */                  
   char*    opts;                   /* all options from commandline */

   /* --- PIPELINE OPTIONS --- */
   int      pipeline_mode;          /* which workflow pipeline to use */
   int      search_mode;            /* alignment search mode */
   int      verbose_level;          /* levels of verbosity */
   char*    tmp_folderpath;         /* location to build a temporary work folder */
   bool     tmp_remove;             /* should temp folder be removed at the end? */
   bool     filter_on;              /* filter thresholds enforced, or let all through? */
   int      is_compo_bias;          /* should composition bias filter be applied? */

   /* --- SEARCH/RANGE OPTIONS --- */
   /* type of searches */
   int      qt_search_space;        /* which queries vs which targets? */
   /* for specified range of targets/queries in file */
   RANGE    t_range;                /* start-end range of targets in file (inclusive) */
   RANGE    q_range;                /* start-end range of queries in file (inclusive) */
   RANGE    list_range;             /* start-end range of hitlist */

   /* --- DEBUG OPTIONS --- */
   bool     is_debug;               /* determines whether debug statements appear */
   char*    dbg_folderpath;         /* location for debugging */
   bool     enforce_warnings;       /* if error is caught, force close? */
   bool     adjust_mmseqs_alns;     /* if mmseqs alignments are out-of-bounds, should we 

   /* --- INPUT --- */
   /* file paths */
   char*    t_filepath;             /* filepath to target (hmm or fasta) file */
   char*    q_filepath;             /* filepath to query (fasta) file */
   /* index paths */
   bool     is_indexpath;           /* is an index file supplied? */
   char*    t_indexpath;            /* index filepath for quick access of target (hmm) file */
   char*    q_indexpath;            /* index filepath for quick access of query (fasta) file */
    /* target/query metadata */
   int      t_filetype;             /* enumerated FILETYPE of target file */
   int      q_filetype;             /* enumerated FILETYPE of query file */
   /* mmseqs-plus search (input) */
   char*    mmseqs_res_filepath;    /* filepath to mmseqs .m8 results file */
   /* simple hitlist (input) */
   char*    hitlist_filepath;       /* filepath to simple hitlist */
   /* database size */
   int      t_dbsize;               /* number of targets in database */
   int      q_dbsize;               /* number of queries in database */

   /* --- OUTPUT --- */
   /* standard output path */
   bool     is_redirect_stdout;     /* are we redirecting stdout? */
   char*    output_filepath;        /* filepath to output results to; "!stdout" => stdout */
    /* standard error path */
   bool     is_redirect_stderr;     /* are we redirecting stdout? */
   char*    error_filepath;         /* filepath to output results to; "!stdout" => stdout */
   /* tblout option/path (modeled after HMMER --tblout) */
   bool     is_tblout;              /* report tblout table? */
   char*    tblout_filepath;        /* filepath to output results; if NULL, doesn't output */
   /* m8out option/path (modeled after MMseqs) */
   bool     is_m8out;               /* report m8out table? */
   char*    m8out_filepath;         /* filepath to output results; if NULL, doesn't output */
   /* myout option/path (my custom output) */
   bool     is_myout;               /* report myout table? */
   char*    myout_filepath;         /* filepath to output results; if NULL, doesn't output */
   /* customized output */
   bool     is_customout;           /* report myout table? */
   char*    customout_filepath;     /* filepath to output results; if NULL, doesn't output */
   bool     custom_fields[15];      /* boolean list of which fields should be reported */ 
  
   /* --- MMSEQS --- */
   int      mmseqs_kmer;               /* kmer length */
   int      mmseqs_split;              /* database split size */
   int      mmseqs_prefilter;          /* double-kmer prefilter kscore */
   int      mmseqs_ungapped_vit;       /* ungapped viterbi */
   float    mmseqs_evalue;             /* e-value gapped viterbi */
   float    mmseqs_pvalue;             /* p-value gapped viterbi */

   /* --- MMORE / FB-PRUNER --- */
   /* if viterbi is precomputed, gives starting and ending coords (single result) */
   COORDS   beg;                    /* beginning coordinates of viterbi alignment */
   COORDS   end;                    /* ending coordinates of viterbi alignment */
   /* cloud search tuning parameters */
   float    alpha;                  /* cloud search: x-drop pruning ratio */
   float    beta;                   /* cloud search: x-drop maximum drop before termination */
   int      gamma;                  /* cloud search: number of antidiag passes before pruning  */
   float    mmore_evalue;           /* e-value mmore / fb-pruner */
   float    mmore_pvalue;           /* p-value mmore / fb-pruner */

   /* --- VITERBI / FWDBACK --- */
   float    threshold_vit;          /* threshold for viterbi score */
   float    threshold_fwd;          /* threshold for forward score */
   float    threshold_mmore;        /* threshold for fb-pruner score */
} ARGS;

/* scores */
typedef struct {
   float    nat_sc;  /* in NATS */
   float    pre_sc;  /* in BITS */
   float    seq_sc;  /* in BITS */
   /* bias correction */
   float    bias_correction;
   float    null_sc;
   float    seq_bias;
   float    filter_sc;
   /* final p-value and e-value */
   float    ln_pval;
   float    pval;
   float    eval;
} SCORES;

typedef struct {
   float    program_start;
   float    program_runtime;
   /* indexes */
   float    load_target_index;
   float    load_query_index;
   /* load query/target */
   float    load_target;
   float    load_query;
   /* naive algs */
   float    naive_cloud;         /* naive cloud search */
   /* quadratic algs */
   float    quad_fwd;            /* forward */
   float    quad_bck;            /* backward */
   float    quad_vit;            /* viterbi */
   float    quad_trace;          /* traceback of viterbi */
   float    quad_cloud_fwd;      /* */
   float    quad_cloud_bck;      /* */
   float    quad_merge;          /* */
   float    quad_reorient;       /* */
   float    quad_bound_fwd;      /* bounded forward (requires cloud) */
   float    quad_bound_bck;      /* bounded backward (requires cloud) */
   float    quad_fbpruner_total; /* */
   /* linear algs */   
   float    lin_fwd;             /* forward-backward */
   float    lin_bck;             /* backward */
   float    lin_vit;             /* viterbi */
   float    lin_trace;           /* traceback of viterbi */
   float    lin_cloud_fwd;       /* */
   float    lin_cloud_bck;       /* */
   float    lin_merge;           /* */
   float    lin_reorient;        /* */
   float    lin_bound_fwd;       /* bounded forwarded (requires cloud) */
   float    lin_bound_bck;       /* bounded backward (requires cloud) */
   float    lin_fbpruner_total;  /* */
   /* sparse algs */
   float    sp_fwd;              /* forward-backward */
   float    sp_bck;              /* backward */
   float    sp_vit;              /* viterbi */
   float    sp_trace;            /* traceback of viterbi */
   float    sp_cloud_fwd;        /* */
   float    sp_cloud_bck;        /* */
   float    sp_merge;            /* */
   float    sp_bound_fwd;        /* bounded forwarded (requires cloud) */
   float    sp_bound_bck;        /* bounded backward (requires cloud) */
   float    sp_fbpruner_total;   /* */
} TIMES;

/* hmm location within file */
typedef struct {
   int         id;            /* id number, determined by order in file */
   char*       name;          /* Name of HMM/FASTA in file */
   long        offset;        /* Positional offset of HMM/FASTA into file */
   int         mmseqs_id;     /* id number, referencing mmseqs lookup file */
} F_INDEX_NODE;

/* index for offset locations into a file, searchable by name or id */
typedef struct {
   /* */
   int            N;             /* Number of location index nodes used */
   int            Nalloc;        /* Number of location index nodes allocated */
   F_INDEX_NODE*  nodes;         /* List of nodes for location of each HMM/FASTA in file */    
   /* */
   char*          index_path;    /* Filepath of index file (NULL if built and not loaded) */
   char*          lookup_path;   /* Filepath to mmseqs lookup file (NULL if not used) */
   char*          source_path;   /* Filepath of file being indexed */
   char*          delim;         /* one of more delimiter of header fields */
   /* */
   int            filetype;      /* Type of file being indexed (HMM, FASTA, etc) */
   int            sort_type;     /* Whether the index nodes list has been sorted, and by which field */
   int            mmseqs_names;  /* Whether index is using names from mmseqs lookup */
} F_INDEX;


/* flags for command line arguments */
typedef struct {
   /* */
   char*    name;          /* name of flag */
   /* */
   int      num_args;      /* number of arguments */
   int      data_type;     /* data type of arguments */
   void*    arg_loc;       /* pointer to the location in ARGS to store option argument */
   /* */
   char*    long_flag;     /* long "--" flag */
   char*    short_flag;    /* single character "-" flag */
   char*    desc;          /* description of flag */
} FLAG_CMD;

/* results */
typedef struct {

} M8_RESULTS;

/* results fields */
typedef struct {
   /* result unique id */
   int         result_id;
   /* cloud search lookup index id */
   int         target_id;
   int         query_id;
   /* mmseqs lookup index id */
   int         target_mid;
   int         query_mid;
   /* name */
   char*       target_name;
   char*       query_name;
   /* characteristics  */
   float       perc_id;
   int         aln_len;
   int         mismatch;
   int         gap_openings;
   /* alignment window */
   COORDS      beg;
   COORDS      end;
   /* alignment window */
   int         query_start;
   int         query_end;
   int         target_start;
   int         target_end;
   /* number of cells computed */
   int         cpu_cloud_cells;     /* number of times a cell is computed */
   int         cloud_cells;         /* number of cells in cloud search matrix */
   int         cpu_total_cells;     /* number of times a cell is computed */
   int         total_cells;         /* total number of cells in full search matrix */
   /* alignment */
   ALIGNMENT*  aln;
   char*       cigar;
   /* times */
   double      total;
   TIMES       times;
   /* m8 fields */ 
   float       e_value;
   int         bit_score;
   float       cloud_fwd_sc; 
   /* scores */
   float       vit_natsc;
   float       fwd_natsc;
   float       bck_natsc;
   float       fbpruner_fwd_natsc;
   float       fbpruner_bck_natsc;
   SCORES      final_scores;
   /* time spent getting result */
   float       time;
} RESULT;

/* results and stats */
typedef struct {
   /* results */
   size_t   N;                      /* number of current results in queue */
   size_t   Nalloc;                 /* allocated space in queue */
   RESULT*  data;                   /* result queue */
   /* aggregate stats */
   int      num_searches;           /* total number of searches */
   int      num_hits;               /* total number of searches to pass threshold */
   float    reporting_threshold;    /* eval threshold for input mmseqs */
   /* output */
   int      max_in_queue;           /* number of results to keep in memory before dumping to file */
   char*    filepath;               /* path to write file */
   FILE*    fp;                     /* file to write to */
} RESULTS;

/* substitution/scoring matrix */
typedef struct {
   char*    filename;               /* */
   char*    alph;                   /* */
   int      alph_len;               /* */
   int      map[256];               /* */
   float*   scores;                 /* */
} SCORE_MATRIX;

/* data struct for passing debugger data */
typedef struct {
   /* debug options */
   bool           is_debugging;     /* debug flag determines whether testing debug output is generated */
   int            verbose_level;    /* level of verbosity: amount to output to user */  
   bool           is_embed;         /* whether to embed linear space matrices in quadratic matrices */
   bool           is_viz;           /* whether to fill data viz matrix */
   /* output paths */
   char*          dbgout_dir;       /* folder */
   char*          dbgout_path;      /* filepath */
   FILE*          dbgout_fp;        /* file pointer */
   /* data structs */
   MATRIX_3D*     test_MX;          /* quadratic testing matrix */
   MATRIX_3D*     test_MX3;         /* lin-space testing matrix */
   MATRIX_2D*     cloud_MX;         /* cloud matrix for visualizing */
   MATRIX_2D*     cloud_MX3;        /* cloud matrix for visualizing */
   EDGEBOUNDS*    test_edg;         /* edgebounds for comparing diags to rows */
} DEBUG_KIT;

/* bools of all tasks to be executed by WORKER */
/* flags set by pipeline and then passed to a generic workflow */
typedef struct {
   /* output heatmaps */
   bool     heatmaps;            /* output heatmaps */
   /* mmseqs lookup */
   bool     mmseqs_lookup;       /* if we have a mmseqs name lookup table */
   /* naive algs */
   bool     naive;            
   bool     naive_cloud;         /* naive cloud search */
   /* quadratic algs */
   bool     quadratic;           /* are we running any quadratic-space algorithms? */
   bool     quad_fwd;            /* forward */
   bool     quad_bck;            /* backward */
   bool     quad_vit;            /* viterbi */
   bool     quad_trace;          /* traceback of viterbi */
   bool     quad_bound_fwd;      /* bounded forward (requires cloud) */
   bool     quad_bound_bck;      /* bounded backward (requires cloud) */
   bool     quad_bias_corr;      /* bias correction */
   /* linear algs */
   bool     linear;              /* are we running any linear-space algorithms? */
   bool     lin_fwd;             /* forward-backward */
   bool     lin_bck;             /* backward */
   bool     lin_vit;             /* viterbi */
   bool     lin_trace;           /* traceback of viterbi */
   bool     lin_cloud_fwd;       /* forward cloud search */
   bool     lin_cloud_bck;       /* forward cloud search */
   bool     lin_bound_fwd;       /* bounded forwarded (requires cloud) */
   bool     lin_bound_bck;       /* bounded backward (requires cloud) */
   /* sparse algs */
   bool     sparse;              /* are we running any linear-space algorithms? */
   bool     sparse_fwd;          /* forward-backward */
   bool     sparse_bck;          /* backward */
   bool     sparse_vit;          /* viterbi */
   bool     sparse_trace;        /* traceback of viterbi */
   bool     sparse_bound_fwd;    /* bounded forwarded (requires cloud) */ 
   bool     sparse_bound_bck;    /* bounded backward (requires cloud) */
   bool     sparse_bias_corr;    /* bias correction */
} TASKS;

typedef struct {
   /* naive algs */
   float     naive_bound_fwd;     /* bound forward */
   float     naive_bound_bck;     /* bound backward */
   /* quadratic algs */
   float     quad_fwd;            /* forward */
   float     quad_bck;            /* backward */
   float     quad_vit;            /* viterbi */
   float     quad_trace;          /* traceback of viterbi */
   float     quad_bound_fwd;      /* bound forward */
   float     quad_bound_bck;      /* bound backward */
   /* linear algs */   
   float     lin_fwd;             /* forward-backward */
   float     lin_bck;             /* backward */
   float     lin_vit;             /* viterbi */
   float     lin_trace;           /* traceback of viterbi */
   float     lin_bound_fwd;       /* bound forward */
   float     lin_bound_bck;       /* bound backward */
   /* sparse algs */
   float     sparse_fwd;          /* forward-backward */
   float     sparse_bck;          /* backward */
   float     sparse_vit;          /* viterbi */
   float     sparse_trace;        /* traceback of viterbi */
   float     sparse_bound_fwd;    /* bound forward */ 
   float     sparse_bound_bck;    /* bound backward */
} NAT_SCORES;

/* collection of all tuning parameters for cloud search */
typedef struct {
   float       alpha;         /* x-drop for local,  */
   float       beta;          /* x-drop for global, determines when to terminate search. looser (larger) than alpha. */
   int         gamma;         /* number of traversed antidiags before pruning begins */
} CLOUD_PARAMS;

/* aggregate stats */
typedef struct {
   /* times */
   float    time_start;       /* start of program */
   float    time_end;         /* end of program */
   float    total_time;       /*  */
   /* database sizes */
   int      n_query_db;       /* number of queries in database */
   int      n_target_db;      /* number of targets in database */
   int      n_query_search;   /* number of queries in search list */
   int      n_target_search;  /* number of targets in search list */
   int      nodes_total;      /* number of hmm model nodes */
   int      resides_total;    /* total number of sequence residues */
   int      n_searches;       /* number of individual searches */
   /* threshold passes */
   int      n_passed_prefilter;    /* number passed prefilter (mmseqs only) */
   int      n_passed_vit;          /* number passed viterbi filter */
   int      n_passed_fwd;          /* number passed forward filter */
   int      n_passed_fbpruner;     /* number passed fbpruner filter */
   int      n_passed_queries;      /* number of queries which reported over threshold */
} STATS;

/* TODO: pass list of target/query ids for searching */
typedef struct {
   int            N;       /* list length */
   int            Nalloc;  /* allocate list length  */
   int*           q_ids;   /* query id list */
   int*           t_ids;   /* target id list */
} HITLIST;

/* single domain data */
typedef struct {
   /* envelope range */
   int            env_beg; 
   int            env_end;
   /* alignement range */
   int            aln_beg;
   int            aln_end;
   /* scores */
   float          env_sc;           /* forward score of envelope */
   float          dom_corr;         /* domain correction -> null2 score when calculating per-domain score (in NATS) */ 
   float          dom_bias;         /* domain bias -> logsum(0, log(bg->omega) + dom_corr ) */ 
   float          optacc_sc;        /* optimal accuraccy score: (units: expected # residues correctly aligned) */
   float          bit_sc;           /* overall score in blocks */
   double         lnP;              /* log(p-value) of the bitscore */
   /* thresholds for reporting */
   bool           is_reported;      /* if domain meets reporting threshold */
   bool           is_included;      /* if domain meets inclusion threshold */
} DOMAIN_X;

/* all domains data */
typedef struct {
   /* domains */
   int            N;             /* number of domains used */
   int            Nalloc;        /* number of domains allocated for */
   DOMAIN_X*      domains;       /* domain array */
   /* vector data */
   VECTOR_FLT*    b_tot;         /* cumulative number of times expected to BEGIN at or before q_0 */
   VECTOR_FLT*    e_tot;         /* cumulative number of times expected to END at or before q_0 */
   VECTOR_FLT*    m_occ;         /* probability X emitted by core model */
   VECTOR_FLT*    null2_exp;     /* null2 expectation by character */
   VECTOR_FLT*    null2_sc;      /* null2 scores by position */
   /* working space for computing null2 score */
   MATRIX_2D*     st_freq;       /* normal state frequencies */
   VECTOR_FLT*    sp_freq;       /* special state frequencies */
   VECTOR_FLT*    st_num;        /* normal state counts (for sparse) */
   /* stats */
   float          n_expected;    /* posterior expectation for number of domains */
   int            n_regions;     /* number of regions */
   int            n_domains;     /* number of domains */
   int            n_envelopes;   /* number of envelopes */
   int            n_clustered;   /* number of clusters of domains */
   int            n_overlaps;    /* number of overlaps of domains */
   /* scores */
   float          dom_fwd;       /* domain forward score */
   float          seq_bias;      /* computed using sum of null score */
   
   /* domain reporting thresholds */
   float          rt1;           /* default =  */
   float          rt2;           /* default =  */
   float          rt3;           /* default =  */
} DOMAIN_DEF;

/* TODO: for multi-threading (stored in WORKER object) */
typedef struct {
    /* currently loaded mmseqs or hitlist id */
   int                  mmseqs_id;
   int                  hitlist_id;
   /* current query and target data */
   SEQUENCE*            q_seq;
   SEQUENCE*            t_seq;
   HMM_PROFILE*         t_prof;
   HMM_BG*              hmm_bg;
   /* edgebounds for cloud search */
   EDGEBOUNDS*          edg_fwd;
   EDGEBOUNDS*          edg_bck;
   EDGEBOUNDS*          edg_diag;
   EDGEBOUNDS*          edg_row;
   /* edgebound row object; helper for reorientation of edgebounds */
   EDGEBOUND_ROWS*      edg_rows_tmp;
   /* int vector for cloud search */
   VECTOR_INT*          lb_vec[3];     /* left-bound lookback vector for bounds */
   VECTOR_INT*          rb_vec[3];     /* right-bound lookback vector for bounds */
   /* cloud pruning parameters */
   CLOUD_PARAMS         cloud_params;  /* */
   /* alignment traceback for viterbi */
   ALIGNMENT*           traceback;     /* traceback */
   /* alignment traceback for maximum posterior */
   ALIGNMENT*           trace_post;    /* posterior traceback */
   /* dynamic programming matrices */
   MATRIX_3D*           st_MX;         /* normal state matrix (quadratic space) */
   MATRIX_3D*           st_MX3;        /* normal state matrix (linear space) */
   MATRIX_3D_SPARSE*    st_SMX;        /* normal state matrix (sparse) */
   MATRIX_2D*           sp_MX;         /* special state matrix (quadratic space) */
   MATRIX_3D*           st_cloud_MX;   /* matrix for naive cloud search */
   /* times for tasks */
   TIMES*               times;
   TIMES*               times_raw;
   /* results to output */
   RESULTS*             results; 
   /* current result */
   RESULT*              result;
} WORKER_THREAD;

/* worker contains the necessary data structures to conduct search */
typedef struct {
   /* arguments set by pipeline defaults and set by user */
   ARGS*                args;
   /* boolean list of pipeline tasks to be performed */
   TASKS*               tasks;
   /* output file pointer */
   FILE*                output_fp;
   FILE*                tblout_fp;
   FILE*                m8out_fp;
   FILE*                myout_fp;
   /* indexes of query and target data files */
   F_INDEX*             q_index;
   F_INDEX*             t_index;
   /* database size for converting p-value to e-value */
   int                  q_size;
   int                  t_size;
   /* file pointers to query and target data file */
   FILE*                q_file;
   FILE*                t_file;
   /* results from mmseqs */
   RESULTS*             results_in;
   /* if given a list of query/target pairs */
   HITLIST*             hitlist_in;    /* */
   /* number of searches */
   int                  num_searches;  /* number of searches */
   /* id of currently loaded query and target */
   int                  q_id;          /* currently loaded query id */
   int                  t_id;          /* currently loaded target id */
   /* currently loaded mmseqs or hitlist id */
   int                  mmseqs_id;     /* */
   int                  hitlist_id;    /* */
   /* current query and target data */
   SEQUENCE*            q_seq;         /* */
   SEQUENCE*            t_seq;         /* */
   HMM_PROFILE*         t_prof;        /* */
   HMM_BG*              hmm_bg;        /* */
   /* edgebounds for cloud search */
   EDGEBOUNDS*          edg_fwd;       /* */
   EDGEBOUNDS*          edg_bck;       /* */
   EDGEBOUNDS*          edg_diag;      /* */
   EDGEBOUNDS*          edg_row;       /* */
   /* edgebound row object; helper for reorientation of edgebounds */
   EDGEBOUND_ROWS*      edg_rows_tmp;  /* */
   /* int vector for cloud search */
   VECTOR_INT*          lb_vec[3];     /* */
   VECTOR_INT*          rb_vec[3];     /* */
   /* cloud pruning parameters */
   CLOUD_PARAMS         cloud_params;  /* */
   /* alignment traceback for viterbi */
   ALIGNMENT*           traceback;     /* */
   ALIGNMENT*           trace_post;    /* */
   /* dynamic programming matrices */
   /* quadratic space */
   MATRIX_3D*           st_MX;         /* normal state matrix (quadratic space) */
   MATRIX_3D*           st_MX_fwd;     /* normal state matrix (quadratic space) */
   MATRIX_3D*           st_MX_bck;     /* normal state matrix (quadratic space) */
   MATRIX_3D*           st_MX_post;    /* normal state matrix (quadratic space) */
   /* linear space */
   MATRIX_3D*           st_MX3;        /* normal state matrix (linear space) */
   MATRIX_3D*           st_MX3_fwd;    /* normal state matrix (linear space) */
   MATRIX_3D*           st_MX3_bck;    /* normal state matrix (linear space) */
   /* sparse */
   MATRIX_3D_SPARSE*    st_SMX;        /* normal state matrix (sparse) */
   MATRIX_3D_SPARSE*    st_SMX_fwd;    /* normal states matrix (sparse), exclusive for forward */
   MATRIX_3D_SPARSE*    st_SMX_bck;    /* normal states matrix (sparse), exclusive for backward */
   MATRIX_3D_SPARSE*    st_SMX_post;   /* normal states matrix (sparse), exclusive for backward */
   /* special state */
   MATRIX_2D*           sp_MX;         /* special state matrix */
   MATRIX_2D*           sp_MX_fwd;     /* special state matrix, exclusive for forward */
   MATRIX_2D*           sp_MX_bck;     /* special state matrix, exclusive for backward */
   MATRIX_2D*           sp_MX_post;    /* special state matrix, exclusive for posterior */
   /* testing */
   MATRIX_3D*           st_cloud_MX;   /* matrix for naive cloud search */
   /* posterior data */
   DOMAIN_DEF*          dom_def;       /* domain boundary data */
   /* times for tasks */
   TIMES*               times;         /* current result section runtimes */
   TIMES*               times_raw;     /* */
   TIMES*               times_totals;  /* cumulative section runtimes */
   double               runtime;       /* total program runtime */
   /* raw scores */
   NAT_SCORES*          scores;        /*  */
   /* results to output */
   RESULTS*             results; 
   /* current result */
   RESULT*              result;
   /* clok for taking times */
   CLOCK*               clok;
   /* TODO: thread objects for multi-threading */
   int                  N_threads;
   int                  Nalloc_threads;
   WORKER_THREAD*       threads;
} WORKER;


/* === GLOBAL VARIABLES === */
/* */
extern char*      PIPELINE_NAMES[];
extern void       (*PIPELINES[])(WORKER*);
extern int        PIPELINE_NUM_ARGS[];
extern char*      MODE_NAMES[];
extern char*      VERBOSITY_NAMES[];
extern char*      ALPHABET_NAMES[];
extern int        ALPHABET_LENGTHS[];
extern char*      STATE_NAMES[];
extern char*      STATE_FULL_NAMES[];
extern char*      STATE_CHARS[];
/* input file types and extensions */
extern char*      FILE_TYPE_EXTS[];
extern char*      FILE_TYPE_NAMES[];
extern int        FILE_TYPE_MAP[];
/* alphabetically-ordered amino lookup and reverse lookup tables */
extern char       ALPH_AMINO_CHARS[];
extern char       AA[];
extern int        AA_REV[];
/* background frequencies of null model, normal and log space */
extern double     BG_MODEL[];
extern double     BG_MODEL_log[];
/* commandline arg objects */
extern ARGS*      args;
extern int        num_flag_cmds;
extern FLAG_CMD   COMMAND_OPTS[];
extern char*      DATATYPE_NAMES[];
/* scoring matrix for converting sequences to hmm */
extern SCORE_MATRIX*    bld;
/* debugging data */
extern DEBUG_KIT*       debugger;

/* script locations */
extern char*       MMSEQS_PLUS_SCRIPT;
extern char*       MMSEQS_PLUS_EASY_SCRIPT;
extern char*       FASTA_TO_HMM_SCRIPT;
/* other tool binary locations */
extern char*       MMSEQS_BIN;
extern char*       HMMBUILD_BIN;

#endif /* _STRUCTS_H */
