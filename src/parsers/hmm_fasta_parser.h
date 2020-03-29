/*******************************************************************************
 *  FILE:      hmm_fasta_parser.h
 *  PURPOSE:   Converts a .fasta to HMM_PROFILE.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/


/* ==== STRUCTS === */
/* imported from HMMER, easel */

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



/* === FUNCTIONS === */

/* Converts .fasta file to build a HMM_PROFILE object */
HMM_PROFILE* HMM_PROFILE_Fasta_Parser( char*   _filename_,
                                       long    offset )
{

}