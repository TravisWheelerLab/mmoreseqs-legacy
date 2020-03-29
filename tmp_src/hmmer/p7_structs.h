/*******************************************************************************
 *  FILE:      p7_structs.c
 *  PURPOSE:   Datatypes imported from HMMER.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _P7_STRUCTS_H
#define _P7_STRUCTS_H

/* === CONSTANTS & ENUMS === */

/* Bit flags used in <hmm->flags>: optional annotation in an HMM
 * 
 * Flags marked with ! may not be changed nor used for other meanings,
 * because they're codes used by HMMER2 (and earlier) that must be
 * preserved for reverse compatibility with old HMMER files.
 * 
 * Why use flags? (So I don't ask this question of myself again:)
 *   1. The way we allocate an HMM, we need to know if we're allocating
 *      M-width annotation fields (RF, CS, CA, MAP) before we read the
 *      annotation from the file.
 *   2. Historically, H2 used flags, so we still need to read H2 flags
 *      for backwards compatibility; so we may as well keep using them.
 */
#define p7H_HASBITS (1<<0)    /* obsolete (was: model has log-odds scores)       !*/
#define p7H_DESC    (1<<1)    /* description exists (legacy; xref SRE:J5/114)    !*/
#define p7H_RF      (1<<2)    /* #RF annotation available                        !*/
#define p7H_CS      (1<<3)    /* #CS annotation available                        !*/
#define p7H_XRAY    (1<<4)    /* obsolete (was: structural data available)       !*/
#define p7H_HASPROB (1<<5)    /* obsolete (was: model in probability form)       !*/
#define p7H_HASDNA  (1<<6)    /* obsolete (was: protein HMM->DNA seq params set) !*/
#define p7H_STATS   (1<<7)    /* model has E-value statistics calibrated         !*/
#define p7H_MAP     (1<<8)    /* alignment map is available                      !*/
#define p7H_ACC     (1<<9)    /* accession is available (legacy; xref SRE:J5/114)!*/
#define p7H_GA      (1<<10)   /* gathering thresholds available                  !*/
#define p7H_TC      (1<<11)   /* trusted cutoffs available                       !*/
#define p7H_NC      (1<<12)   /* noise cutoffs available                         !*/
#define p7H_CA      (1<<13)   /* surface accessibilities available               !*/
#define p7H_COMPO   (1<<14)   /* model-specific residue composition available     */
#define p7H_CHKSUM  (1<<15)   /* model has an alignment checksum                  */
#define p7H_CONS    (1<<16)   /* consensus residue line available                 */
#define p7H_MMASK   (1<<17)   /* #MM annotation available                        !*/

#define p7_MAXABET    20      /* maximum size of alphabet (4 or 20)              */

#define p7_NEVPARAM 6   /* number of statistical parameters stored in models                      */
#define p7_NCUTOFFS 6   /* number of Pfam score cutoffs stored in models                          */
#define p7_NOFFSETS 3   /* number of disk offsets stored in models for hmmscan's fast model input */

enum p7_archchoice_e { p7_ARCH_FAST = 0, p7_ARCH_HAND = 1 };
enum p7_wgtchoice_e  { p7_WGT_NONE  = 0, p7_WGT_GIVEN = 1, p7_WGT_GSC    = 2, p7_WGT_PB       = 3, p7_WGT_BLOSUM = 4 };
enum p7_effnchoice_e { p7_EFFN_NONE = 0, p7_EFFN_SET  = 1, p7_EFFN_CLUST = 2, p7_EFFN_ENTROPY = 3, p7_EFFN_ENTROPY_EXP = 4 };

#define p7O_NXSTATES  4    /* special states stored: ENJC                       */
#define p7O_NXTRANS   2         /* special states all have 2 transitions: move, loop */
#define p7O_NTRANS    8    /* 7 core transitions + BMk entry                    */

enum p7o_xstates_e      { p7O_E    = 0, p7O_N    = 1,  p7O_J  = 2,  p7O_C  = 3 };
enum p7o_xtransitions_e { p7O_MOVE = 0, p7O_LOOP = 1 };
enum p7o_tsc_e          { p7O_BM   = 0, p7O_MM   = 1,  p7O_IM = 2,  p7O_DM = 3, p7O_MD   = 4, p7O_MI   = 5,  p7O_II = 6,  p7O_DD = 7 };

/* Relative entropy target defaults:
 * For proteins, hmmbuild's effective sequence number calculation
 * aims to achieve a certain relative entropy per match emission.
 * (= average score per match emission).
 * These are empirically tuned constants,
 */
#define p7_ETARGET_AMINO  0.59 /* bits,  from the work of Steve Johnson. */
#define p7_ETARGET_DNA    0.62 /* bits,  from the work of Travis Wheeler and Robert Hubley. */
#define p7_ETARGET_OTHER  1.0  /* bits */ /* if you define your own alphabet, set this */

/* Indices for special state types in the length model, gm->xsc[x][]
 */
enum p7p_xstates_e { 
  p7P_E = 0,
  p7P_N = 1,
  p7P_J = 2,
  p7P_C = 3
};
#define p7P_NXSTATES 4

/* Indices of Plan7 main model state transitions, hmm->t[k][] */
enum p7h_transitions_e {
  p7H_MM = 0,
  p7H_MI = 1,
  p7H_MD = 2,
  p7H_IM = 3,
  p7H_II = 4,
  p7H_DM = 5,
  p7H_DD = 6 
};
#define p7H_NTRANSITIONS 7

/* Indices for transitions from the length modeling scores gm->xsc[][x]
 */
enum p7p_xtransitions_e {
  p7P_LOOP = 0,
  p7P_MOVE = 1
};
#define p7P_NXTRANS 2

/* === STRUCTS === */

/* Structure: MIXDCHLET
 * 
 * A mixture Dirichlet density, usually used as a prior 
 * for a multinomial model (turning count vectors into probability
 * parameters).
 */
typedef struct {
 /*::cexcerpt::dirichlet_mixdchlet::begin::*/
  double  *pq;       /* mixture coefficients pq[0..N-1]          */
  double **alpha;               /* Dirichlet params alpha[0..N-1][0..K-1]   */
  int      N;        /* number of mixtures, e.g. 9 for Sjolander */
  int      K;        /* alphabet size, e.g. 20                   */
 /*::cexcerpt::dirichlet_mixdchlet::end::*/
} ESL_MIXDCHLET;

/* */
typedef struct p7_prior_s {
  ESL_MIXDCHLET *tm;    /*  match transitions */
  ESL_MIXDCHLET *ti;    /* insert transitions */
  ESL_MIXDCHLET *td;    /* delete transitions */
  ESL_MIXDCHLET *em;    /*  match emissions   */
  ESL_MIXDCHLET *ei;    /* insert emissions   */
} P7_PRIOR;

/* */
typedef struct p7_hmm_s {
  /*::cexcerpt::plan7_core::begin::*/
  int     M;                    /* length of the model (# nodes)                           */
  float **t;                    /* transition prob's. t[(0),1..M][0..p7H_NTRANSITIONS-1]   */
  float **mat;                  /* match emissions.  mat[1..M][0..K-1]                     */ 
  float **ins;                  /* insert emissions. ins[1..M][0..K-1]                     */
  /*::cexcerpt::plan7_core::end::*/

  /* Annotation. Everything but <name> is optional. Flags are set when
   * optional values are set. All the char *'s are proper nul-terminated
   * strings, not just arrays. (hmm->map is an int array).
   */
  char    *name;                 /* name of the model                     (mandatory)      */ /* String, \0-terminated   */
  char    *acc;                    /* accession number of model (Pfam)      (optional: NULL) */ /* String, \0-terminated   */
  char    *desc;                 /* brief (1-line) description of model   (optional: NULL) */ /* String, \0-terminated   */
  char    *rf;                   /* reference line from alignment 1..M    (p7H_RF)         */ /* String; 0=' ', M+1='\0' */
  char    *mm;                   /* model mask line from alignment 1..M   (p7H_MM)         */ /* String; 0=' ', M+1='\0' */
  char    *consensus;               /* consensus residue line        1..M    (p7H_CONS)       */ /* String; 0=' ', M+1='\0' */
  char    *cs;                   /* consensus structure line      1..M    (p7H_CS)         */ /* String; 0=' ', M+1='\0' */
  char    *ca;                  /* consensus accessibility line  1..M    (p7H_CA)         */ /* String; 0=' ', M+1='\0' */

  char    *comlog;               /* command line(s) that built model      (optional: NULL) */ /* String, \0-terminated   */
  int      nseq;           /* number of training sequences          (optional: -1)   */
  float    eff_nseq;             /* effective number of seqs (<= nseq)    (optional: -1)   */
  int      max_length;           /* upper bound length, all but 1e-7 prob (optional: -1)   */
  char    *ctime;          /* creation date                         (optional: NULL) */
  int     *map;                    /* map of alignment cols onto model 1..M (p7H_MAP)        */ /* Array; map[0]=0 */
  uint32_t checksum;             /* checksum of training sequences        (p7H_CHKSUM)     */
  float    evparam[p7_NEVPARAM]; /* E-value params                        (p7H_STATS)      */
  float    cutoff[p7_NCUTOFFS];  /* Pfam score cutoffs                    (p7H_{GA,TC,NC}) */
  float    compo[p7_MAXABET];    /* model bg residue comp                 (p7H_COMPO)      */

  off_t    offset;               /* HMM record offset on disk                              */
  const ESL_ALPHABET *abc;       /* ptr to alphabet info (hmm->abc->K is alphabet size)    */
  int      flags;                /* status flags                                           */
} P7_HMM;

/*  */
typedef struct p7_bg_s {
  float   *f;     /* null1 background residue frequencies [0..K-1]: set at initialization    */
  float    p1;    /* null1's transition prob: p7_bg_SetLength() sets this from target seq L  */

  ESL_HMM *fhmm;  /* bias filter: p7_bg_SetFilter() sets this, from model's mean composition */

  float    omega; /* the "prior" on null2/null3: set at initialization (one omega for both null types)  */

  const ESL_ALPHABET *abc; /* reference to alphabet in use: set at initialization             */
} P7_BG;

/* */
typedef struct p7_profile_s {
  float  *tsc;          /* transitions  [0.1..M-1][0..p7P_NTRANS-1], hand-indexed  */
  float **rsc;          /* emissions [0..Kp-1][0.1..M][p7P_NR], hand-indexed       */
                        /*           [ AMINO ][ POS  ][ MATCH OR INSERT ]          */
  float   xsc[p7P_NXSTATES][p7P_NXTRANS]; /* special transitions [NECJ][LOOP,MOVE] */

  int     mode;         /* configured algorithm mode (e.g. p7_LOCAL)               */ 
  int     L;      /* current configured target seq length                    */
  int     allocM; /* max # of nodes allocated in this structure              */
  int     M;      /* number of nodes in the model                            */
  int     max_length;   /* calculated upper bound on emitted seq length            */
  float   nj;     /* expected # of uses of J; precalculated from loop config */

  /* Info, most of which is a copy from parent HMM:                                       */
  char  *name;       /* unique name of model                                   */
  char  *acc;        /* unique accession of model, or NULL                     */
  char  *desc;                  /* brief (1-line) description of model, or NULL           */
  char  *rf;                    /* reference line from alignment 1..M; *rf=0 means unused */
  char  *mm;                    /* modelmask line           1..M; *ref=0: unused     */
  char  *cs;                    /* consensus structure line      1..M, *cs=0 means unused */
  char  *consensus;     /* consensus residues to display in alignments, 1..M      */
  float  evparam[p7_NEVPARAM];   /* parameters for determining E-values, or UNSET          */
  float  cutoff[p7_NCUTOFFS];    /* per-seq/per-domain bit score cutoffs, or UNSET         */
  float  compo[p7_MAXABET];   /* per-model HMM filter composition, or UNSET             */

  /* Disk offset information for hmmpfam's fast model retrieval                           */
  off_t  offs[p7_NOFFSETS];     /* p7_{MFP}OFFSET, or -1                                  */

  off_t  roff;                  /* record offset (start of record); -1 if none            */
  off_t  eoff;                  /* offset to last byte of record; -1 if unknown           */

  const ESL_ALPHABET *abc; /* copy of pointer to appropriate alphabet                */
} P7_PROFILE;

/* temp */
typedef struct p7_oprofile_s {
  /* MSVFilter uses scaled, biased uchars: 16x unsigned byte vectors                 */
  __m128i **rbv;         /* match scores [x][q]: rm, rm[0] are allocated      */
  __m128i **sbv;         /* match scores for ssvfilter                        */
  uint8_t   tbm_b;    /* constant B->Mk cost:    scaled log 2/M(M+1)       */
  uint8_t   tec_b;    /* constant E->C  cost:    scaled log 0.5            */
  uint8_t   tjb_b;    /* constant NCJ move cost: scaled log 3/(L+3)        */
  float     scale_b;    /* typically 3 / log2: scores scale to 1/3 bits      */
  uint8_t   base_b;            /* typically +190: offset of uchar scores            */
  uint8_t   bias_b;    /* positive bias to emission scores, make them >=0   */

  /* ViterbiFilter uses scaled swords: 8x signed 16-bit integer vectors              */
  __m128i **rwv;    /* [x][q]: rw, rw[0] are allocated  [Kp][Q8]         */
  __m128i  *twv;    /* transition score blocks          [8*Q8]           */
  int16_t   xw[p7O_NXSTATES][p7O_NXTRANS]; /* NECJ state transition costs            */
  float     scale_w;            /* score units: typically 500 / log(2), 1/500 bits   */
  int16_t   base_w;             /* offset of sword scores: typically +12000          */
  int16_t   ddbound_w;    /* threshold precalculated for lazy DD evaluation    */
  float     ncj_roundoff;  /* missing precision on NN,CC,JJ after rounding      */

  /* Forward, Backward use IEEE754 single-precision floats: 4x vectors               */
  __m128 **rfv;         /* [x][q]:  rf, rf[0] are allocated [Kp][Q4]         */
  __m128  *tfv;          /* transition probability blocks    [8*Q4]           */
  float    xf[p7O_NXSTATES][p7O_NXTRANS]; /* NECJ transition costs                   */

  /* Our actual vector mallocs, before we align the memory                           */
  __m128i  *rbv_mem;
  __m128i  *sbv_mem;
  __m128i  *rwv_mem;
  __m128i  *twv_mem;
  __m128   *tfv_mem;
  __m128   *rfv_mem;
  
  /* Disk offset information for hmmpfam's fast model retrieval                      */
  off_t  offs[p7_NOFFSETS];     /* p7_{MFP}OFFSET, or -1                             */

  /* Disk offset bookkeeping for h3f:                                                */
  off_t  roff;                  /* record offset (start of record); -1 if none       */
  off_t  eoff;                  /* offset to last byte of record; -1 if unknown      */

  /* Information, annotation copied from parent profile:                             */
  char  *name;      /* unique name of model                              */
  char  *acc;      /* unique accession of model, or NULL                */
  char  *desc;                  /* brief (1-line) description of model, or NULL      */
  char  *rf;                    /* reference line           1..M; *ref=0: unused     */
  char  *mm;                    /* modelmask line           1..M; *ref=0: unused     */
  char  *cs;                    /* consensus structure line 1..M, *cs=0: unused      */
  char  *consensus;    /* consensus residues for ali display, 1..M          */
  float  evparam[p7_NEVPARAM];   /* parameters for determining E-values, or UNSET     */
  float  cutoff[p7_NCUTOFFS];   /* per-seq/per-dom bit cutoffs, or UNSET             */
  float  compo[p7_MAXABET];  /* per-model HMM filter composition, or UNSET        */
  const ESL_ALPHABET *abc;  /* copy of ptr to alphabet information               */

  /* Information about current configuration, size, allocation                       */
  int    L;      /* current configured target seq length              */
  int    M;      /* model length                                      */
  int    max_length;    /* upper bound on emitted sequence length            */
  int    allocM;    /* maximum model length currently allocated for      */
  int    allocQ4;    /* p7_NQF(allocM): alloc size for tf, rf             */
  int    allocQ8;    /* p7_NQW(allocM): alloc size for tw, rw             */
  int    allocQ16;    /* p7_NQB(allocM): alloc size for rb                 */
  int    mode;      /* currently must be p7_LOCAL                        */
  float  nj;      /* expected # of J's: 0 or 1, uni vs. multihit       */

  int    clone;                 /* this optimized profile structure is just a copy   */
                                /* of another profile structre.  all pointers of     */
                                /* this structure should not be freed.               */
} P7_OPROFILE;

/* */
typedef struct p7_trace_s {
  int    N;    /* length of traceback                       */
  int    nalloc;        /* allocated length of traceback             */
  char  *st;      /* state type code                   [0..N-1]*/
  int   *k;    /* node index; 1..M if M,D,I; else 0 [0..N-1]*/
  int   *i;    /* pos emitted in dsq, 1..L; else 0  [0..N-1]*/
  float *pp;      /* posterior prob of x_i; else 0     [0..N-1]*/
  int    M;    /* model length M (maximum k)                */
  int    L;    /* sequence length L (maximum i)             */

  /* The following section is data generated by "indexing" a trace's domains */
  int   ndom;     /* number of domains in trace (= # of B or E states) */
  int  *tfrom,   *tto;  /* locations of B/E states in trace (0..tr->N-1)     */
  int  *sqfrom,  *sqto; /* first/last M-emitted residue on sequence (1..L)   */
  int  *hmmfrom, *hmmto;/* first/last M state on model (1..M)                */
  int   ndomalloc;   /* current allocated size of these stacks            */
} P7_TRACE;

/* */
typedef struct p7_builder_s {
  /* Model architecture                                                                            */
  enum p7_archchoice_e arch_strategy;    /* choice of model architecture determination algorithm   */
  float                symfrac;           /* residue occ thresh for fast architecture determination */
  float                fragthresh;   /* if L <= fragthresh*alen, seq is called a fragment      */

  /* Relative sequence weights                                                                     */
  enum p7_wgtchoice_e  wgt_strategy;     /* choice of relative sequence weighting algorithm        */
  double               wid;       /* %id threshold for BLOSUM relative weighting            */

  /* Effective sequence number                                                                     */
  enum p7_effnchoice_e effn_strategy;    /* choice of effective seq # determination algorithm      */
  double               re_target;    /* rel entropy target for effn eweighting, if set; or -1.0*/
  double               esigma;       /* min total rel ent parameter for effn entropy weights   */
  double               eid;       /* %id threshold for effn clustering                      */
  double               eset;      /* effective sequence number, if --eset; or -1.0          */

  /* Run-to-run variation due to random number generation                                          */
  ESL_RANDOMNESS      *r;           /* RNG for E-value calibration simulations                */
  int                  do_reseeding;    /* TRUE to reseed, making results reproducible            */

  /* E-value parameter calibration                                                                 */
  int                  EmL;                /* length of sequences generated for MSV fitting          */
  int                  EmN;            /* # of sequences generated for MSV fitting               */
  int                  EvL;                /* length of sequences generated for Viterbi fitting      */
  int                  EvN;            /* # of sequences generated for Viterbi fitting           */
  int                  EfL;            /* length of sequences generated for Forward fitting      */
  int                  EfN;            /* # of sequences generated for Forward fitting           */
  double               Eft;            /* tail mass used for Forward fitting                     */

  /* Choice of prior                                                                               */
  P7_PRIOR            *prior;          /* choice of prior when parameterizing from counts        */
  int                  max_insert_len;

  /* Optional: information used for parameterizing single sequence queries                         */
  ESL_SCOREMATRIX     *S;      /* residue score matrix                                   */
  ESL_DMATRIX         *Q;           /* Q->mx[a][b] = P(b|a) residue probabilities             */
  double               popen;           /* gap open probability                                   */
  double               pextend;          /* gap extend probability                                 */

  double               w_beta;    /*beta value used to compute W (window length)   */
  int                  w_len;     /*W (window length)  explicitly set */

  const ESL_ALPHABET  *abc;       /* COPY of alphabet                                       */
  char errbuf[eslERRBUFSIZE];            /* informative message on model construction failure      */
} P7_BUILDER;


#endif /* _P7_STRUCTS_H */