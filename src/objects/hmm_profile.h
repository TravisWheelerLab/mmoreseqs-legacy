/*******************************************************************************
 *  @file hmm_profile.c
 *  @brief HMM_PROFILE Object
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

#ifndef _HMM_PROFILE_H
#define _HMM_PROFILE_H

// /* === OBJECTS === */
// typedef enum {
//    AMINO, 
//    DNA
// } ALPHABET;

// typedef struct {
//    float param1;
//    float param2;
// } DIST_PARAM;

// typedef struct {
//    /* NORMAL STATE PROBABILITIES */
//    float          match[NUM_AMINO];          /* match emission probabilities for each amino acid */
//    float          insert[NUM_AMINO];         /* insert emission probabilities for each amino acid */
//    /* [0]A  [1]C  [2]D  [3]E  [4]F  [5]G  [6]H  [7]I  [8]K  [9]L
//       [10]M [11]N [12]P [13]Q [14]R [15]S [16]T [17]V [18]W [19]Y */
//    float          trans[NUM_TRANS_STATES];   /* transition state probabilities (default same as COMPO) */
//    /* [0]m->m [1]m->i [2]m->d [3]i->m [4]i->i [5]d->m [6]d->d [7]b->m */
// } HMM_NODE;

// typedef struct {
//    float          freq[NUM_AMINO];           /* hard-coded background residue frequencies for each amino acid */
//    float          compo[NUM_AMINO];          /* background residue frequencies of the given hmm model */
//    float          insert[NUM_AMINO];         /* insert emission probabilities for each amino acid (uniform across positions) */
//     [0]A  [1]C  [2]D  [3]E  [4]F  [5]G  [6]H  [7]I  [8]K  [9]L
//       [10]M [11]N [12]P [13]Q [14]R [15]S [16]T [17]V [18]W [19]Y 
//    float          trans[NUM_TRANS_STATES];   /* transition state probabilities (default same as COMPO) */
//    /* [0]m->m [1]m->i [2]m->d [3]i->m [4]i->i [5]d->m [6]d->d [7]b->m */

//    /* SPECIAL STATE PROBABILITIES */
//    float          spec[NUM_SPECIAL_STATES][NUM_SPECIAL_TRANS];
//    /* [0]N [1]E [2]C [3]J */
//    /* [0]LOOP  [1]MOVE */
//    int            num_J;
// } HMM_BG;

// typedef struct {
//    int            N;                   /* profile length (number of nodes) */
//    int            alph_leng;           /* alphabet length: AMINO = 20, DNA = 4 */
//    /* profile settings */
//    int            isLocal; 
//    int            isMultihit;     
//    /* */
//    char           *filename;
//    /* meta data */
//    char           *name; 
//    char           *acc; 
//    char           *desc;
//    char           *alph;               /* alphabet type () */;      
//    /* distribution parameters for scoring */
//    DIST_PARAM     *msv_dist; 
//    DIST_PARAM     *viterbi_dist; 
//    DIST_PARAM     *forward_dist; 
//    /* models */
//    HMM_BG         *bg_model;           /* background composition */
//    HMM_NODE       *hmm_model;          /* array of position specific probabilities */
// } HMM_PROFILE;


/* === FUNCTIONS === */
/* Constructor */
HMM_PROFILE* HMM_PROFILE_Create();
/* Destructor */
void HMM_PROFILE_Destroy( HMM_PROFILE* prof );

/* Set Textfield to HMM_PROFILE field */
void HMM_PROFILE_Set_TextField( char** prof_field, 
                                char*  text );
/* Set HMM Model Length and allocate memory for nodes */
void HMM_PROFILE_Set_Model_Length( HMM_PROFILE* prof, 
                                   int          length );
/* Set alphabet (DNA or AMINO ACID) for HMM_PROFILE */
void HMM_PROFILE_Set_Alphabet( HMM_PROFILE* prof, 
                               char*        alph_name );
/* Set Distribution Parameters for HMM_PROFILE */
void HMM_PROFILE_Set_Distribution_Params( HMM_PROFILE* prof, 
                                          float        param1, 
                                          float        param2, 
                                          char*        dist_name );

/* Output HMM_PROFILE to FILE POINTER */
void HMM_PROFILE_Dump( HMM_PROFILE* prof, 
                       FILE*        fp );

#endif /* _HMM_PROFILE_H */
