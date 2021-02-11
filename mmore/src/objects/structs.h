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

/* === MACROS === */
#include "structs_macros.h"
#include "structs_consts.h"
#include "structs_funcs.h"
#include "structs_buildopts.h"

/* === DECLARE HIDDEN FUNCTIONS (for older compilers) === */
/* (*** c99 standard hides this function from <string.h> ***) */
// extern char* strdup(const char*);
// extern ssize_t getline(char **lineptr, size_t *n, FILE *stream);

/* === EASEL === */
#include "easel.h"
#include "esl_hmm.h"
#include "esl_alphabet.h"
#include "esl_sq.h"
#include "esl_stopwatch.h"

/* === STRUCTS === */

/* --- BASIC DATA TYPES --- */
/* template datatype */
#define XXX          int
/* alias primitive datatypes */
#define BOOL         bool
#define INT          int
#define FLT          float
#define CHAR         char 
#define DBL          double
#define STR          char*
/* datatype for matrices */
#define DATA         float

/* tool for describing data */
typedef struct {
   DATATYPE    type;
   STR         name;
   STR         desc;
   size_t      size;
} DATA_DESC;

/* generic datatype that can hold most basic datatypes */
typedef union {
   int      i;
   float    f;
   double   d;
   long     l;
   char*    s;
   char     c;
   bool     b;
   void*    p;
} GENERIC;

/* generic data with specified type */
typedef struct {
   GENERIC        data;
   DATATYPE       type;
} GEN;

/* coordinates in matrix */
typedef struct {
   int      q_0;       /* row index */
   int      t_0;       /* col index */
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

/* string object  */
typedef struct {
   int      N;
   char*    data;
} XSTR;

/* given cell of alignment */
typedef struct {
   int      q_0;       /* index in query */
   int      t_0;       /* index in target */
   int      st;        /* state at index */
} TRACE;

/* === VECTORS === */

/* default size of vectors when created */
#define VECTOR_INIT_SIZE 128

/* vector of template data */
typedef struct {
   XXX*        data;     /* array of data type */
   size_t      N;        /* current length of array in use */
   size_t      Nalloc;   /* current length of array allocated */
}  VECTOR_XXX;

/* vector of characters */
typedef struct {
   CHAR*    data;             /* array of data type */
   size_t   N;                /* current length of array in use */
   size_t   Nalloc;           /* current length of array allocated */
}  VECTOR_CHAR;

/* vector of strings */
typedef struct {
   STR*     data;             /* array of data type */
   size_t   N;                /* current length of array in use */
   size_t   Nalloc;           /* current length of array allocated */
}  VECTOR_STR;

/* vector of integers */
typedef struct {
   INT*     data;             /* array of data type */
   size_t   N;                /* current length of array in use */
   size_t   Nalloc;           /* current length of array allocated */
}  VECTOR_INT;

/* vector of floats */
typedef struct {
   FLT*     data;             /* array of data type */
   size_t   N;                /* current length of array in use */
   size_t   Nalloc;           /* current length of array allocated */
}  VECTOR_FLT;

/* vector of doubles */
typedef struct {
   DBL*     data;             /* array of data type */
   size_t   N;                /* current length of array in use */
   size_t   Nalloc;           /* current length of array allocated */
}  VECTOR_DBL;

/* vector of edgebounds  */
typedef struct {
   BOUND*   data;    /* array of data type */
   size_t   N;       /* current length of array in use */
   size_t   Nalloc;  /* current length of array allocated */
}  VECTOR_BOUND;

/* vector of range structs */
typedef struct {
   RANGE*   data;     /* array of data type */
   size_t   N;        /* current length of array in use */
   size_t   Nalloc;   /* current length of array allocated */
}  VECTOR_RANGE;

/* vector of trace structs */
typedef struct {
   TRACE*   data;     /* array of data type */
   size_t   N;        /* current length of array in use */
   size_t   Nalloc;   /* current length of array allocated */
}  VECTOR_TRACE;

/* vector of data structs */
typedef struct {
   DATA*    data;     /* array of data type */
   size_t   N;        /* current length of array in use */
   size_t   Nalloc;   /* current length of array allocated */
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

/* === INPUT/OUTPUT === */

/* file manager for opening/closing files */
typedef struct {
   FILE*             fp;            /* pointer to file to be open/closed */
   STR               filename;      /* path to file */
   STR               mode;          /* which mode is file to be opened in? (e.g. "a", "r", "w", "b", "r+") */
   bool              is_open;       /* is file currently opened? */
   bool              is_eof;        /* is file pointer at the end of file? */
   long int          prv_pos;       /* position of start of current line in file */
   long int          cur_pos;       /* position of beginning of next line in file */
} FILER;

/* string buffer for reading/writing to file */
typedef struct {    
   /* current line buffer */
   VECTOR_CHAR*      line;          /* last line read to buffer */
   STR               buffer_ptr;    /* pointer into <line> buffer */
   int               line_len;      /* number of characters in current line; -1 if eof */
   /* split line buffer by words/fields */
   VECTOR_CHAR*      split_line;    /* copy of line for splitting */ 
   STR               field;         /* pointer to next field in <split_line> */
   VECTOR_INT*       field_offsets; /* offsets to the head of all fields in <split_line> */
   STR               delim;         /* field separator / delimiters */
   int               field_idx;     /* current index in <field_offsets> */
} BUFFER;

/* reader for parsing files */
typedef struct {
   /* file data */
   FILE*             fp;            /* file to be read */
   char*             filename;      /* string to file location */
   long int          file_size;     /* size of file (in chars) */
   bool              is_open;       /* whether file is open or not */
   bool              is_eof;        /* is reader at end of file? */
   size_t            line_count;    /* line number in file */
   size_t            lines_read;    /* number of lines read */
   /* buffer data */
   BUFFER*           buffer;        /* buffer for data */
   /* position in file */
   long int          pos_prv;       /* position of start of previous line */
   long int          pos_cur;       /* position of start of current line */
   long int          pos_nxt;       /* position of start of next line */
} READER;

/* reader for parsing files */
typedef struct {
   /* file data */
   FILER*            file;          /* file to be written to */
   BUFFER*           buffer;        /* buffer for data */
   /* meta data */
   size_t            num_writes;    /* number of writes to file */
} WRITER;

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
   RANGE          Q_range;    /* subrange of query that rows should be allocated for */
   int            Q;          /* size of query */
   int            T;          /* size of target */
   /* data */
   BOUND*         rows;       /* array of bounds, each for a specific row */
} EDGEBOUND_ROWS;

/* alignment for viterbi traceback */
typedef struct {
   /* dimensions of embedded matrix */
   int            Q;             /* length of query sequence */
   int            T;             /* length of target hmm profile */

   /* full model traceback from T->S */
   int            full_len;      /* length of complete T->S alignment */
   VECTOR_TRACE*  traces;        /* list of all (state,i,j) TRACES in ALIGNMENT */
   VECTOR_FLT*    scores;        /* list of scores at each position in the traceback */
   /* discrete core model alignments (separated by jumps) */
   int            num_alns;      /* number of distinct alignment */
   VECTOR_INT*    tr_beg;        /* index of every alignment begin (B state) in traceback */
   VECTOR_INT*    tr_end;        /* index of every alignment end (E state) in traceback */
   VECTOR_FLT*    tr_score;      /* cumulative score for every corresponding <tr_beg,tr_end> alignment in traceback */

   /* best alignment */
   int            best_idx;      /* index of best alignment */
   size_t         aln_len;       /* alignment length */
   /* trace endpoints */
   int            beg;           /* current beginning index in traces */
   int            end;           /* current end index in traces */
   /* counts of hits, mismatches, and gaps */
   int            num_gaps;      /* number of gaps in alignment */
   int            num_misses;    /* number of misses in alignment */
   int            num_matches;   /* number of matches */
   float          perc_id;       /* percent identity */

   /* string representations of alignment */
   /* mmseqs' cigar-style alignments */
   bool           is_cigar_aln;  /* has cigar-style alignment been constructed? */
   VECTOR_CHAR*   cigar_aln;     /* ex: 379M173M2D41M2D65M6I21 (see MMseqs userguide pg.35) */
   /* hmmer-style alignments */
   bool           is_hmmer_aln;  /* has hmmer-style alignment been constructed? */
   VECTOR_CHAR*   target_aln;    /* |qdrfLePqcrilDlevWdqCYfRWlPvLeikgGGqpq| */
   VECTOR_CHAR*   center_aln;    /* |++  L    ++  l+vW qCYfRW+P  e+ +G    | */
   VECTOR_CHAR*   query_aln;     /* |DEDLLSMDATVKALSVWTQCYFRWVPKAEVVNGDPGT| */
   VECTOR_CHAR*   state_aln;     /* |MMMMIIIDIDMM...                      | */
} ALIGNMENT;

/* clock for timing events (wrapper for squid stopwatch above) */
typedef struct {
   /* current time */
   double         start;            /* captures start time */
   double         stop;             /* captures stop time */
   double         duration;         /* captures difference between start and stop */
   /* previous times */
   double         program_start;    /* time captured when clock is created */
   VECTOR_DBL*    stamps;
   #if ( CLOCK_TYPE == CLOCK_EASEL ) 
   /* wrapper for easel watch */
   ESL_STOPWATCH* esl_stopwatch;
   double         esl_total;
   #endif
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
   VECTOR_CHAR*   consensus;        /* consensus sequence */
   int            is_consensus;     /* has consensus been built */
   /* main model */
   HMM_COMPO*     bg_model;         /* background composition */
   HMM_NODE*      hmm_model;        /* array of position specific probabilities */
   /* submodel */
   int            N_full;           /* profile length of full model */
   HMM_NODE*      hmm_model_full;   /* array of position at start of full model */
   /* database */
   int            N_tprofs;         /* number of target hmm profile models */
   int            N_qseqs;          /* number of query sequences */
} HMM_PROFILE;

/* Sequence */
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
   int            q_0;           /* current query position iterator */
   int            t_0;           /* current target position iterator */
   RANGE          r_0;           /* current edgebound iterator for retrieving next row range */
} MATRIX_3D_SPARSE;

/* dynamic programming matrix for computing algs */
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
   DPMX_MODE            mode;          /* type of matrix currently being stored */
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
   bool     tmp_remove;             /* should temp files/folders be removed at the end? */

   /* --- SEARCH/RANGE OPTIONS --- */
   /* type of searches */
   int      qt_search_space;        /* which queries vs which targets? */
   /* for specified range of targets/queries in file */
   RANGE    t_range;                /* start-end range of targets in file (inclusive) */
   RANGE    q_range;                /* start-end range of queries in file (inclusive) */
   RANGE    list_range;             /* start-end range of hitlist */

   /* --- DEBUG OPTIONS --- */
   bool     is_use_local_tools;     /* whether to system installed tools or local project tools */
   bool     is_recycle_mx;          /* whether to recycle <fwd> and <bck> matrices for computing <post> and <optacc> */
   bool     is_debug;               /* determines whether debug statements appear */
   char*    dbg_folderpath;         /* location for debugging */
   bool     enforce_warnings;       /* if error is caught, force close? */
   bool     adjust_mmseqs_alns;     /* if mmseqs alignments are out-of-bounds, should we truncate alignment? */

   /* --- INPUT --- */
   /* file paths */
   char*    t_filepath;             /* filepath to target (hmm, msa, or fasta) file */
   char*    q_filepath;             /* filepath to query (fasta) file */
   char*    t_mmseqs_filepath;      /* filepath to mmseqs target (hhm, msa, mm_msa or fasta) file */
   /* target/query metadata */
   bool           is_guess_filetype;      /* whether to use guessing tool to find file type */
   FILE_TYPE      t_filetype;             /* enumerated FILETYPE of target file */
   FILE_TYPE      q_filetype;             /* enumerated FILETYPE of query file */
   FILE_TYPE      t_mmseqs_filetype;      /* enumerated FILETYPE of mmseqs target file */
   /* index paths */
   bool     is_indexpath;           /* is an index file supplied? */
   char*    t_indexpath;            /* index filepath for quick access of target (hmm) file */
   char*    q_indexpath;            /* index filepath for quick access of query (fasta) file */
   /* database size */
   int      t_dbsize;               /* number of targets in database */
   int      q_dbsize;               /* number of queries in database */

   /* --- INTERRIM OUTPUT --- */
   /* mmseqs-plus search (input) */
   char*    mmseqs_m8_filepath;    /* filepath to mmseqs .m8 results file */
   /* simple hitlist (input) */
   char*    hitlist_filepath;       /* filepath to simple hitlist */
   /* optional output */
   bool     is_mmseqs_p2sout;       /* output profile-to-sequence mmseqs results? */
   char*    mmseqs_p2s_filepath;    /* filepath to output results; if NULL, doesn't output */
   bool     is_mmseqs_s2sout;       /* output profile-to-sequence mmseqs results? */
   char*    mmseqs_s2s_filepath;    /* filepath to output results; if NULL, doesn't output */

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
   /* mydomout option/path (my custom output for domains) */
   bool     is_mydomout;            /* report mydomout table? */
   char*    mydomout_filepath;      /* filepath to output results; if NULL, doesn't output */
   /* mytimeout option/path (my custom output for times) */
   bool     is_mytimeout;            /* report mytimeout table? */
   char*    mytimeout_filepath;      /* filepath to output results; if NULL, doesn't output */
   /* mythreshout option/path (my custom output for times) */
   bool     is_mythreshout;            /* report mythreshout table? */
   char*    mythreshout_filepath;      /* filepath to output results; if NULL, doesn't output */
   
   /* customized output */
   bool     is_customout;           /* report myout table? */
   char*    customout_filepath;     /* filepath to output results; if NULL, doesn't output */
   bool     custom_fields[15];      /* boolean list of which fields should be reported */ 
  
   /* --- TASK OPTIONS --- */
   bool     is_run_pruned;          /* run pruned forward backward */
   bool     is_run_full;            /* run full forward backward */
   bool     is_run_domains;         /* run domain search */
   bool     is_compo_bias;          /* should composition bias filter be applied? */

   /* --- MMSEQS --- */
   int      mmseqs_hits_per_search;    /* maximum number of alignments allowed to be reported per target/query search */ 
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

   /* --- THRESHOLDS --- */
   bool     is_run_filter;          /* filter thresholds on mmore enforced, or let all through? */
   float    threshold_pre;          /* threshold for prefilter score */
   float    threshold_vit;          /* threshold for viterbi score */
   float    threshold_fwd;          /* threshold for forward score */
   float    threshold_cloud;        /* threshold for cloud search score */
   float    threshold_bound_fwd;    /* threshold for bound forward score */
   float    threshold_mmore;        /* threshold for mmore score */
} ARGS;

/* scores */
typedef struct {
   /* --- alignment scores --- */
   /* nat-scores */
   float    prefilter_natsc;        /* mmseqs double-kmer prefilter natscore */
   float    viterbi_natsc;          /* viterbi natscore */
   float    cloud_natsc;            /* cloud search composite natscore */
   float    bound_fwdback_natsc;    /* bound forward backward natscore */
   float    fwdback_natsc;          /* full forward backward natscore */
   /* bit-scores */
   float    prefilter_bitsc;        /* mmseqs double-kmer prefilter bitscore */
   float    viterbi_bitsc;          /* viterbi bitscore */
   float    cloud_bitsc;            /* cloud search composite bitscore */
   float    bound_fwdback_bitsc;    /* bound forward backward bitscore */
   float    fwdback_bitsc;          /* full forward backward bitscore */
   /* e-values */ 
   float    prefilter_eval;         /* mmseqs double-kmer prefilter eval (current DNE) */
   float    viterbi_eval;           /* viterbi eval */
   float    cloud_eval;             /* cloud search composite eval */
   float    bound_fwdback_eval;     /* bound forward backward eval */
   float    fwdback_eval;           /* full forward backward eval */

   /* --- bias correction scores --- */
   /* nat-scores */
   float    null_omega_natsc;       /* prior probability of no bias natscore */
   float    null1_hmm_bias_natsc;   /* null1 hmm model bias natscore */ 
   float    null2_seq_bias_natsc;   /* null2 sequence bias natscore */
   float    filtersc_natsc;         /* null1 filter score (current not used) */
   /* bit-scores */
   float    null_omega_bitsc;       /* prior probability of no bias bitscore */
   float    null1_hmm_bias_bitsc;   /* null1 hmm model bias bitscore */ 
   float    null2_seq_bias_bitsc;   /* null2 sequence bias bitscore */
   float    filtersc_bitsc;         /* null1 filter score (current not used) */

   /* final scores */
   float    nat_sc;     /* final forward score w/o correction (in NATS) */
   float    pre_sc;     /* final forward w/ null1 bias correction (in BITS) */
   float    seq_sc;     /* final forward w/ null1 and null2 bias correction (in BITS) */
   float    sum_sc;     /* final forward constructed by summing over all domain seq_sc (in BITS) */
   float    pval;       /* final forward pval (best of seq_sc and sum_sc) */
   float    eval;       /* final forward eval (best of seq_sc and sum_sc) */
} SCORES;

/* times for subroutines during program lifetime */
typedef struct {
   /* --- overall --- */
   float    program_start;       /* start time of program */
   float    program_end;         /* end time of program */
   float    program;             /* runtime of the program */

   /* --- pre-pipeline loop tasks --- */
   float    parse_args;          /* parse args */
   float    load_target_index;   /* load or build target index */
   float    load_query_index;    /* load or build query index */

   /* --- pipeline loop overall --- */
   float    loop_start;          /* start time of loop iteration */
   float    loop_end;            /* end time of loop iteration */
   float    loop;                /* runtime of loop iteration */

   /* --- pipeline loop tasks --- */
   /* load target or query */
   float    load_target;         /* load next target hmm profile for loop */
   float    load_query;          /* load next query seqeunce for loop */
   /* naive algs */
   float    naive_cloud;         /* naive cloud search */
   /* quadratic algs */
   float    quad_fwd;            /* forward */
   float    quad_bck;            /* backward */
   float    quad_vit;            /* viterbi */
   float    quad_trace;          /* viterbi traceback */
   float    quad_cloud_fwd;      /* cloud search forward */
   float    quad_cloud_bck;      /* cloud search backward */
   float    quad_merge;          /* edgebounds merge/union */
   float    quad_reorient;       /* edgebounds reorient */
   float    quad_bound_fwd;      /* bound forward */
   float    quad_bound_bck;      /* bound backward */
   float    quad_posterior;      /* posterior computations */
   float    quad_optacc;         /* optimal accuracy and traceback */
   /* linear algs */   
   float    lin_fwd;             /* forward */
   float    lin_bck;             /* backward */
   float    lin_vit;             /* viterbi */
   float    lin_trace;           /* traceback of viterbi */
   float    lin_cloud_fwd;       /* cloud search forward */
   float    lin_cloud_bck;       /* cloud search backward */
   float    lin_merge;           /* edgebounds merge */
   float    lin_reorient;        /* edgebounds reorient */
   float    lin_bound_fwd;       /* bound forward */
   float    lin_bound_bck;       /* bound backward */
   /* sparse algs */
   float    sp_build_mx;         /* build sparse matrix */
   float    sp_fwd;              /* forward-backward */
   float    sp_bck;              /* backward */
   float    sp_vit;              /* viterbi */
   float    sp_trace;            /* traceback of viterbi */
   float    sp_cloud_fwd;        /* cloud search forward */
   float    sp_cloud_bck;        /* cloud search backward */
   float    sp_bound_fwd;        /* bound forward */
   float    sp_bound_bck;        /* bound backward */
   float    sp_posterior;        /* posterior computations */
   float    sp_decodedom;        /* decode domains */
   float    sp_biascorr;         /* null2 bias correction */
   float    sp_optacc;           /* optimal accuracy and traceback */

   /* --- domain loop --- */
   /* domain-specific (overall times for all domains for given search) */
   float    dom_start;           /* start dom time */
   float    dom_end;             /* end dom time */
   float    dom_total;           /* total dom runtime */
   float    dom_bound_fwd;       /* bound forward */
   float    dom_bound_bck;       /* bound backward */
   float    dom_posterior;       /* domain posterior */
   float    dom_biascorr;        /* null2 bias correction */
   float    dom_optacc;          /* optimal accuracy and traceback */
} TIMES;

/* all scores for all types of algs */
typedef struct {
   /* naive algs */
   float    naive_bound_fwd;        /* bound forward */
   float    naive_bound_bck;        /* bound backward */
   /* quadratic algs */
   float    quad_fwd;               /* forward */
   float    quad_bck;               /* backward */
   float    quad_vit;               /* viterbi */
   float    quad_cloud_fwd;         /* cloud forward search */
   float    quad_cloud_bck;         /* cloud backward search */
   float    quad_bound_fwd;         /* bound forward */
   float    quad_bound_bck;         /* bound backward */
   /* linear algs */   
   float    lin_fwd;                /* forward-backward */
   float    lin_bck;                /* backward */
   float    lin_vit;                /* viterbi */
   float    lin_cloud_fwd;          /* cloud forward search */
   float    lin_cloud_bck;          /* cloud backward search */
   float    lin_bound_fwd;          /* bound forward */
   float    lin_bound_bck;          /* bound backward */
   /* sparse algs */
   float    sparse_fwd;             /* forward-backward */
   float    sparse_bck;             /* backward */
   float    sparse_vit;             /* viterbi */
   float    sparse_bound_fwd;       /* bound forward */ 
   float    sparse_bound_bck;       /* bound backward */

   /* composition bias */
   float    null_omega;             /* presumed prior probability of no bias (currently hardcoded in HMM_BG) */
   float    null1_filtersc;         /* (currently unused) */
   float    null1_hmm_bias;         /* hmm model composition bias */
   float    null2_seq_bias;         /* sequence composition bias */
   /* final scores */
   float    pre_score;              /* score that only accounts in null1 bias */
   float    seq_score;              /* score that accounts for null1 and null2 bias */
   float    sum_score;              /* score that sums score over all domains */
   float    pval;                   /* score converted to p-value */
   float    eval;                   /* score converted to e-value */

   /* threshold scores */
   float    threshold_vit;          /* threshold for viterbi */
   float    threshold_cloud_max;    /* threshold for cloud search max */
   float    threshold_cloud_compo;  /* threshold for cloud search composite */
   float    threshold_bound_max;    /* threshold for bound fwdback max */
   float    threshold_dom_max;      /* threshold for max domain fwdback */
   float    threshold_dom_compo;    /* threshold for composite domain fwdback */
} ALL_SCORES;

/* hmm location within file */
typedef struct {
   int         id;            /* id number, determined by order in file */
   char*       name;          /* Name of HMM/FASTA in file */
   long        offset;        /* Positional offset of HMM/FASTA into file */
   int         mmseqs_id;     /* id number, referencing mmseqs lookup file */
} F_INDEX_NODE;

/* index for offset locations into a file, searchable by name or id */
typedef struct {
   /* node vector */
   int            N;             /* Number of location index nodes used */
   int            Nalloc;        /* Number of location index nodes allocated */
   F_INDEX_NODE*  nodes;         /* List of nodes for location of each HMM/FASTA in file */    
   /* sources */
   char*          index_path;    /* Filepath of index file (NULL if built and not loaded) */
   char*          lookup_path;   /* Filepath to mmseqs lookup file (NULL if not used) */
   char*          source_path;   /* Filepath of file being indexed */
   char*          delim;         /* one of more delimiter of header fields */
   /* meta data */
   int            filetype;      /* Type of file being indexed (HMM, FASTA, etc) */
   int            sort_type;     /* Whether the index nodes list has been sorted, and by which field */
   int            mmseqs_names;  /* Whether index is using names from mmseqs lookup */
} F_INDEX;

/* descriptor for command line arguments */
typedef struct {
   /* output data for user */
   char*    name;          /* name of flag */
   /* program data */
   int      num_args;      /* number of arguments */
   int      data_type;     /* data type of arguments */
   void*    arg_loc;       /* pointer to the location in ARGS to store option argument */
   void*    arg_bool;      /* toggle boolean? */
   /* input data for user */
   char*    long_flag;     /* long "--" flag */
   char*    short_flag;    /* single character "-" flag */
   char*    desc;          /* description of flag (for --help) */
} ARG_OPT;

/* TODO: formatted data */
/* descriptor for data field */
typedef struct {
   STR         name;       /* name of field */
   STR         desc;       /* description of field */
   DATATYPE    type;       /* datatype of field */
   void*       data_loc;   /* data location in worker to retrieve */
} DATA_FIELD;

/* descriptor for data format */
typedef struct {
   STR         name;          /* name of data format */
   STR         desc;          /* description of data format */
   int         num_fields;    /* number of fields in data format */
   DATA_FIELD* field_desc;    /* descriptor for all data fields */
} DATA_FORMAT;

/* formatted results */
typedef struct {
   DATA_FORMAT*   format;           /* description of data format */
   int            total;            /* total number of results processed */
   int            N_in_queue;       /* total number of results in queue */
   int            N_max_in_queue;   /* number of results allowed in queue before outputting */
   DATA_FIELD*    fields;           /* field data */
} FORMATTED_RESULT;

/* m8 result data entry */
typedef struct {
   /* result unique id (for mmore pipeline, this is simply the position in mmseqs output) */
   int         result_id;     /* unique id for this result (generally simple ordering of result in file) */
   /* target/query id */
   int         target_id;     /* target hmm profile id */
   int         query_id;      /* query sequence id */
   /* target/query name */
   char*       target_name;   /* target hmm profile name */
   char*       query_name;    /* query sequence name */
   /* details */
   float       perc_id;       /* percentage of identical matches */
   int         aln_len;       /* alignment length */
   int         mismatch;      /* number of mismatches */
   int         gap_openings;  /* number of gap openings */
   int         q_beg;         /* beginning of query */
   int         q_end;         /* ending of query */
   int         t_beg;         /* beginnning of target */
   int         t_end;         /* ending of query */
   float       eval;          /* e-value / expected value */
   float       bitsc;         /* bit score */
   /* optional field */
   STR         cigar_aln;     /* MMSEQS-style, cigar alignment */
} M8_RESULT;

/* total mmseqs m8 data array */
typedef struct {
   /* file data */
   STR         filename;
   /* data entries */
   int         N;
   int         Nalloc;
   int         N_max;
   int         total;
   M8_RESULT*  data;
   /* stats */
   int         num_searches;
   int         num_hits;
} M8_RESULTS;

/* hitlist entry */
typedef struct {
   /* results id */
   int         result_id;     /* unique id for this result (generally simple ordering of result in file) */
   /* target/query id */
   int         target_id;     /* target hmm profile id */
   int         query_id;      /* query sequence id */      
} HITLIST_RESULT;

/* TODO: pass list of target/query ids for searching */
typedef struct {
   int            N;       /* list length */
   int            Nalloc;  /* allocate list length  */
   int*           q_ids;   /* query id list */
   int*           t_ids;   /* target id list */
} HITLIST_RESULTS;

/* results fields */
typedef struct {
   /* result unique id (for mmore pipeline, this is simply the ordering in mmseqs output) */
   int         result_id;              /* result id */
   int         target_id;              /* target id (not in use) */
   int         query_id;               /* query id (not in use) */
   /* target/query name */
   STR         target_name;            /* target name */
   STR         query_name;             /* query name */
   /* target/query alignment ranges */
   RANGE       target_range;           /* range of target in alignment */
   RANGE       query_range;            /* range of query in alignment */
   /* alignment */
   ALIGNMENT*  vit_trace;              /* viterbi alignment */
   ALIGNMENT*  post_trace;             /* posterior optimal accuracy alignment */
   /* scores */
   ALL_SCORES  scores;                 /* all possible scores produced by algorithms */
   SCORES      final_scores;           /* final scores for reporting */
   /* number of cells computed */
   int         cpu_cloud_cells;        /* number of times a cell in cloud edgebounds */
   int         cloud_cells;            /* number of cells in cloud search matrix */
   int         cpu_total_cells;        /* number of times a cell is computed */
   int         total_cells;            /* total number of cells in full search matrix */
   float       perc_cells;             /* percent of total cells contained in cloud search */
   /* threshold filters (score, cutoff threshold, and boolean check) */
   float       threshold_prefilter;    /* minimum required score to pass prefilter threshold */
   float       score_prefilter;        /* score achieved by prefilter */ 
   bool        is_passed_prefilter;    /* whether prefilter score passed filter */
   float       threshold_viterbi;      /* minimum required score to pass viterbi threshold */
   float       score_viterbi;          /* score achieved by viterbi */ 
   bool        is_passed_viterbi;      /* whether viterbi score passed filter */
   float       threshold_cloud;        /* minimum required score to pass cloud threshold */
   float       score_cloud;            /* score achieved by cloud */ 
   bool        is_passed_cloud;        /* whether cloud score passed filter */
   float       threshold_fwdback;      /* minimum required score to pass fwdback threshold */
   float       score_fwdback;          /* score achieved by fwdback */ 
   bool        is_passed_fwdback;      /* whether fwdback score passed filter */
} RESULT;

/* results and stats */
typedef struct {
   /* results */
   size_t   N;                      /* number of current results in queue */
   size_t   Nalloc;                 /* allocated space in queue */
   RESULT*  data;                   /* result queue */
   int      max_in_queue;           /* number of results to keep in memory before dumping to file */
   /* aggregate data */
   int      num_hits;
   int      num_searches;
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
   bool     naive;               /* are we running any naive algorithms? */
   bool     naive_cloud;         /* naive cloud search */
   /* quadratic algs */
   bool     quadratic;           /* are we running any quadratic-space algorithms? */
   bool     quad_fwd;            /* forward */
   bool     quad_bck;            /* backward */
   bool     quad_vit;            /* viterbi */
   bool     quad_trace;          /* viterbi traceback */
   bool     quad_bound_fwd;      /* bound forward */
   bool     quad_bound_bck;      /* bound backward */
   bool     quad_bias_corr;      /* bias correction */
   /* linear algs */
   bool     linear;              /* are we running any linear-space algorithms? */
   bool     lin_fwd;             /* forward-backward */
   bool     lin_bck;             /* backward */
   bool     lin_vit;             /* viterbi */
   bool     lin_trace;           /* viterbi traceback */
   bool     lin_cloud_fwd;       /* forward cloud search */
   bool     lin_cloud_bck;       /* forward cloud search */
   bool     lin_bound_fwd;       /* bound forward */
   bool     lin_bound_bck;       /* bound backward */
   /* sparse algs */
   bool     sparse;              /* are we running any linear-space algorithms? */
   bool     sparse_fwd;          /* forward-backward */
   bool     sparse_bck;          /* backward */
   bool     sparse_vit;          /* viterbi */
   bool     sparse_trace;        /*  viterbi traceback */
   bool     sparse_bound_fwd;    /* bound forward */ 
   bool     sparse_bound_bck;    /* bound backward */
   bool     sparse_bias_corr;    /* bias correction */
} TASKS;

/* collection of all tuning parameters for cloud search */
typedef struct {
   float       alpha;         /* x-drop for antidiagonal score. */
   float       beta;          /* x-drop for global score, determines when to terminate search. Looser than alpha. */
   int         gamma;         /* number of traversed antidiags before pruning begins */
} CLOUD_PARAMS;

/* aggregate stats */
typedef struct {
   /* database sizes */
   int      n_query_db;             /* number of queries in database */
   int      n_target_db;            /* number of targets in database */
   int      n_query_search;         /* number of queries in search list */
   int      n_target_search;        /* number of targets in search list */
   int      n_nodes_total;          /* total number of hmm model nodes */
   int      n_resides_total;        /* total number of sequence residues */
   int      n_searches;             /* number of individual searches */
   /* threshold passes */
   int      n_passed_prefilter;     /* number passed prefilter (mmseqs only) */
   int      n_passed_viterbi;       /* number passed viterbi filter */
   int      n_passed_cloud;         /* number passed forward filter */
   int      n_passed_fwdback;       /* number passed fbpruner filter */
   /* expected threshold passes */
   int      n_expected_prefilter;   /* number expected to pass prefilter (mmseqs only) */
   int      n_expected_viterbi;     /* number expected to pass viterbi filter */
   int      n_expected_cloud;       /* number expected to pass cloud filter */
   int      n_expected_fwdback;     /* number expected to pass forward filter */
   /* reported */
   int      n_reported_searches;    /* number of searches reported */
   int      n_reported_domains;     /* number of domains reported */
   int      n_reported_targets;     /* number of target profiles reported */
   int      n_reported_queries;     /* number of query sequences reported */
} STATS;

/* single domain data */
typedef struct {
   /* envelope range */
   int            env_beg; 
   int            env_end;
   /* alignment range */
   RANGE          Q_range;
   RANGE          T_range;
   /* alignments */
   ALIGNMENT*     traceback;
   /* scores and e-values */
   float          fwd_sc;
   float          bck_sc;
   float          null2_seq_bias;
   /* final scores */
   float          env_sc;           /* forward score of domain envelope */
   float          dom_corr;         /* domain correction -> null2 score when calculating per-domain score (in NATS) */ 
   float          dom_bias;         /* domain bias -> logsum(0, log(bg->omega) + dom_corr ) */ 
   float          optacc_sc;        /* optimal accuraccy score: (units: expected # residues correctly aligned) */
   float          bit_sc;           /* overall score in bits */
   double         lnP;              /* log(p-value) of the bitscore */
} DOMAIN_X;

/* all domains data */
typedef struct {
   /* domain array */
   int            N;             /* number of domains used */
   int            Nalloc;        /* number of domains allocated */
   DOMAIN_X*      domains;       /* domain array */
   VECTOR_RANGE*  dom_ranges;    /* domain start-end points */
   VECTOR_FLT*    dom_fwdsc;     /* domain forward score */
   VECTOR_FLT*    dom_bias;      /* domain sequence bias correction (via null2) */
   /* current domain */
   int            idx;           /* index of current domain */
   DOMAIN_X       domain;        /* temp for current domain */

   /* top-scoring domain */
   int            best;          /* index of highest scoring domain */
   float          best_sc;       /* corrected score of best domain */
   float          best_presc;    /* pre-score of best domain */
   float          best_fwdsc;    /* forward score of best domain */
   float          best_bias;     /* compo bias of best domain */
   RANGE          best_range;    /* query range of best domain */
   
   /* composite sum-scoring over all domains */
   float          nullsc;           /* null bias */
   float          null_omega;       /* prior probability of no bias */
   float          null1_hmm_bias;   /* null1 background hmm profile bias */
   float          null2_seq_bias;   /* null2 background sequence bias */
   float          dom_sumsc;        /* all domains sumscore */
   float          dom_sumbias;      /* all domains biases */
   
   /* working data for finding domains */
   VECTOR_FLT*    b_tot;         /* cumulative probability of having reached BEGIN state */
   VECTOR_FLT*    e_tot;         /* cumulative probability of having reached END state */
   VECTOR_FLT*    m_occ;         /* cumulative probability of being in the CORE states (M,I,D) */
   /* working data for computing null2 composition bias */
   VECTOR_FLT*    null2_exp;     /* null2 expectation by character */
   VECTOR_FLT*    null2_sc;      /* null2 scores by position */
   /* working data for computing null2 score */
   MATRIX_2D*     st_freq;       /* normal state frequencies */
   VECTOR_FLT*    sp_freq;       /* special state frequencies */
   VECTOR_FLT*    st_num;        /* normal state counts (for sparse) */
   /* working data for finding alignments */
   ALIGNMENT*     align;         /* alignment */
   VECTOR_TRACE*  trace_aln;     /* traceback of alignment */
   VECTOR_CHAR*   cigar_aln;     /* cigar alignment (ex: "TBBBMIDEEEEC") */
   /* working data for temporary domain edgebounds */
   EDGEBOUNDS*    edg;           /* edgebounds for domain */

   /* domain-specific stats */
   int            n_residues;    /* number of residues covered by domains */
   int            n_domains;     /* number of domains */
   int            n_regions;     /* number of regions */
   int            n_envelopes;   /* number of envelopes */
   int            n_clustered;   /* number of clustered regions */
   int            n_overlaps;    /* number of overlapping regions */

   /* domain region threshold cutoffs (for finding domain edges) */
   float          rt1;           /* default region threshold */
   float          rt2;           /* default region threshold */
   float          rt3;           /* default region threshold */
} DOMAIN_DEF;

/* tools for running script scripts */
typedef struct {
   /* main script */
   STR               tool;          /* script interpreter (e.g. "bash", "python", etc) */
   STR               script_path;   /* script filepath */
   /* environmental variables */
   VECTOR_STR*       env_names;     /* environmental variable name */
   VECTOR_STR*       env_values;    /* environmental variable value */
   /* script arguments */
   VECTOR_STR*       arg_flags;     /* script argument flag */
   VECTOR_STR*       arg_values;    /* script argument value */
   /* command executor */
   VECTOR_STR*       command;       /* stores final command to be passed to shell */
} SCRIPTRUNNER;

/* TODO: for multi-threading (stored in WORKER object) */
typedef struct {
   /* --- thread identifier --- */
   int                  thread_id;


} WORKER_THREAD;

/* worker contains the necessary data structures to conduct search */
typedef struct {
   /* --- pipeline --- */
   // PIPELINE*            pipeline;
   /* --- pipeline settings --- */
   ARGS*                args;          /* arguments set by pipeline defaults and set by user */
   TASKS*               tasks;         /* boolean list of pipeline tasks to be performed */
   /* --- pipeline tools --- */
   CLOCK*               timer;         /* stopwatch for timing events */
   /* --- script running tool --- */
   SCRIPTRUNNER*        runner;

   /* --- reader/writer & buffer --- */
   READER*              reader;              /* general use reader */
   // WRITER*              writer;              /* general use writer */
   BUFFER*              buffer;              /* general use read buffer */

   /* --- input files --- */
   FILER*               q_seq_file;          /* query sequence file */
   FILER*               t_prof_file;         /* target hmm profile file */
   FILER*               q_index_file;        /* query index file */
   FILER*               t_index_file;        /* target index file */
   FILER*               mmseqs_file;         /* mmseqs *.m8 format file */
   FILER*               hitlist_file;        /* hitlist csv file */

   /* --- output files --- */
   FILER*               output_file;         /* standard output */
   FILER*               tblout_file;         /* HMMER-style .tblout format output */
   FILER*               m8out_file;          /* MMseqs/BLAST-style .m8 format output */
   FILER*               myout_file;          /* Custom per-search tsv output */
   FILER*               mydomout_file;       /* Custom per-domain tsv output */
   FILER*               mytimeout_file;      /* Runtime summary output */
   FILER*               mythreshout_file;    /* Threshold passage output */

   /* --- input data --- */
   /* m8 results from mmseqs */
   M8_RESULTS*          mmseqs_data;      /* mmseqs .m8 data  */
   HITLIST_RESULTS*     hitlist_data;     /* simple list of target/query pairs (not implemented) */
   /* indexes of query and target data files */
   F_INDEX*             q_index;          /* file index of <q_file> */
   F_INDEX*             t_index;          /* file index of <t_file> */

   /* --- output data --- */
   /* aggregate statistics */
   STATS*               stats;         /* statistics */
   /* results to output (if we are going to output results in batches) */
   RESULTS*             results;       /* results array */
   RESULT*              result;        /* current result */
   /* times for tasks */
   TIMES*               times;         /* current result section runtimes */
   TIMES*               times_totals;  /* cumulative section runtimes */
   /* scores (moved to results) */
   // ALL_SCORES*          scores;        /* scores in NATS */
   // ALL_SCORES*          evals;         /* scores in evals */

   /* --- working data --- */
   /* current query and target data */
   SEQUENCE*            q_seq;         /* query sequence model data */
   SEQUENCE*            t_seq;         /* target sequence model data */
   HMM_PROFILE*         t_prof;        /* target hmm profile model data */
   HMM_BG*              hmm_bg;        /* hmm background model */
   /* edgebounds for cloud search */
   EDGEBOUNDS*          edg_fwd;       /* edgebounds for forward cloud search */
   EDGEBOUNDS*          edg_bck;       /* edgebounds for backward cloud search */
   EDGEBOUNDS*          edg_diag;      /* merged cloud search by antidiagonal */
   EDGEBOUNDS*          edg_row;       /* merged cloud search by row */
   EDGEBOUND_ROWS*      edg_rows_tmp;  /* temporary edgebound row object; helper for reorientating */
   /* int vector for cloud search */
   VECTOR_INT*          lb_vec[3];     /* left bounds for building cloud edgebounds */
   VECTOR_INT*          rb_vec[3];     /* right bounds for building cloud edgebounds */
   /* cloud pruning parameters */
   CLOUD_PARAMS         cloud_params;  /* parameters for cloud search */
   /* alignment traceback for viterbi */
   ALIGNMENT*           trace_vit;     /* traceback for viterbi */
   ALIGNMENT*           trace_post;    /* traceback for posterior */
   /* dynamic programming matrices */
   /* quadratic space matrices */
   MATRIX_3D*           st_MX;         /* normal state matrix (quadratic space) */
   MATRIX_3D*           st_MX_fwd;     /* normal state matrix (quadratic space), exclusive for forward */
   MATRIX_3D*           st_MX_bck;     /* normal state matrix (quadratic space), exclusive for backward */
   MATRIX_3D*           st_MX_post;    /* normal state matrix (quadratic space), exclusive for posterior */
   MATRIX_3D*           st_MX_optacc;  /* normal state matrix (quadratic space), exclusive for optimal accuracy */
   /* linear space matrices */
   MATRIX_3D*           st_MX3;        /* normal state matrix (linear space) */
   MATRIX_3D*           st_MX3_fwd;    /* normal state matrix (linear space), exclusive for forward */
   MATRIX_3D*           st_MX3_bck;    /* normal state matrix (linear space), exclusive for backward */
   /* sparse matrices */
   MATRIX_3D_SPARSE*    st_SMX;        /* normal state matrix (sparse) */
   MATRIX_3D_SPARSE*    st_SMX_fwd;    /* normal states matrix (sparse), exclusive for forward */
   MATRIX_3D_SPARSE*    st_SMX_bck;    /* normal states matrix (sparse), exclusive for backward */
   MATRIX_3D_SPARSE*    st_SMX_post;   /* normal states matrix (sparse), exclusive for posterior */
   MATRIX_3D_SPARSE*    st_SMX_optacc; /* normal states matrix (sparse), exclusive for optimal accuracy */
   /* special state matrices */
   MATRIX_2D*           sp_MX;         /* special state matrix */
   MATRIX_2D*           sp_MX_fwd;     /* special state matrix, exclusive for forward */
   MATRIX_2D*           sp_MX_bck;     /* special state matrix, exclusive for backward */
   MATRIX_2D*           sp_MX_post;    /* special state matrix, exclusive for posterior */
   MATRIX_2D*           sp_MX_optacc;  /* special state matrix, exclusive for optimal accuracy */
   /* posterior data */
   DOMAIN_DEF*          dom_def;       /* domain boundary data */

   /* --- loop data & variables --- */
   /* search id */
   RANGE                search_rng;       /* search range */
   int                  search_id;        /* search id (position in mmseqs source file) */
   int                  n_searches;       /* total number of searches in file */
   int                  n_searches_run;   /* total number of searches run */
   /* mmseqs variables */
   int                  mmseqs_id;        /* current mmseqs id (position in mmseqs loaded data) */
   M8_RESULT*           mmseqs_cur;       /* current mmseqs entry */
   M8_RESULT*           mmseqs_prv;       /* previous mmseqs entry */
   /* hitlist variables */
   int                  hitlist_id;       /* currently loaded hitlist entry from input hitlist results */
   HITLIST_RESULT*      hitlist_cur;      /* current hitlist entry */
   HITLIST_RESULT*      hitlist_prv;      /* previous hitlist entry */
   /* current iteration query/target id */
   int                  q_id;             /* current query id */
   int                  t_id;             /* current target id */
   int                  q_id_prv;         /* previous query id */
   int                  t_id_prv;         /* previous target id */
   /* current iteration query/target name */
   STR                  q_name;           /* current query name */
   STR                  t_name;           /* current target name */
   STR                  q_name_prv;       /* previous query name */
   STR                  t_name_prv;       /* previous target name */
   
   /* --- multi-threading --- */
   int                  N_threads;        /* Number of threads in use */
   int                  Nalloc_threads;   /* Number of threads allocated for */
   WORKER_THREAD*       threads;          /* worker threads array */
} WORKER;

typedef struct {
   size_t               N;
   size_t               Nalloc;
   WORKER*              workers;
} MASTER_WORKER;

/* pipeline descriptors */
typedef struct {
   char*          name;                   /* name of function  */
   STATUS_FLAG    (*func)(WORKER*);       /* pointer to main pipeline function */
   int            num_main_args;          /* number of main args  */
   // STATUS_FLAG    (*set_args)(WORKER*);   /* pointer to set default args function */
   // STATUS_FLAG    (*set_tasks)(WORKER*);  /* pointer to set default tasks function */
} PIPELINE;

/* === GLOBAL VARIABLES === */
/* pipeline */
extern PIPELINE      PIPELINES[];
extern char*         MODE_NAMES[];
extern char*         VERBOSITY_NAMES[];
extern char*         ALPHABET_NAMES[];
extern int           ALPHABET_LENGTHS[];
extern char*         STATE_NAMES[];
extern char*         STATE_FULL_NAMES[];
extern char*         STATE_CHARS[];
/* input file types and extensions */
extern char*         FILE_TYPE_EXTS[];
extern char*         FILE_TYPE_NAMES[];
extern int           FILE_TYPE_MAP[];
/* alphabetically-ordered amino lookup and reverse lookup tables */
extern char          ALPH_AMINO_CHARS[];
extern char          AA[];
extern int           AA_REV[];
/* background frequencies of null model, normal and log space */
extern double        BG_MODEL[];
extern double        BG_MODEL_log[];
/* commandline arg objects */
extern int           num_flag_cmds;
extern ARG_OPT       COMMAND_OPTS[];
extern char*         DATATYPE_NAMES[];
/* scoring matrix for converting sequences to hmm */
extern SCORE_MATRIX* bld;

/* root directory */
extern char*         ROOT_DIR;
/* script locations */
extern char*         MMSEQS_PLUS_SCRIPT;
extern char*         MMSEQS_PLUS_EASY_SCRIPT;
extern char*         FASTA_TO_HMM_SCRIPT;
/* other tool binary locations */
extern char*         MMSEQS_BIN;
extern char*         HMMER_BIN;
extern char*         MMORE_BIN;

/* debugging data */
extern DEBUG_KIT*    debugger;

#endif /* _STRUCTS_H */
