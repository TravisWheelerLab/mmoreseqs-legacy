/*******************************************************************************
 *  - FILE:   alignment.c
 *  - DESC:    ALIGNMENT Object.
 *          Used for storing Viterbi or Forward/Backward Alignments.
 *******************************************************************************/

#ifndef _ALIGNMENT_H
#define _ALIGNMENT_H

/*! FUNCTION:  ALIGNMENT_Create()
 *  SYNOPSIS:  Constructs <aln>, allocates memory and returns pointer.
 */
ALIGNMENT* ALIGNMENT_Create();

/*! FUNCTION:  ALIGNMENT_Destroy()
 *  SYNOPSIS:  Destroys <aln>, frees memory and return NULL pointer.
 */
ALIGNMENT* ALIGNMENT_Destroy(ALIGNMENT* aln);

/*! FUNCTION:  ALIGNMENT_Reuse()
 *  SYNOPSIS:  Wipes <aln>'s old data for reuse, sets dimensions <Q x T>.
 */
void ALIGNMENT_Reuse(ALIGNMENT* aln, int Q, int T);

/*! FUNCTION:  ALIGNMENT_GetSize()
 *  SYNOPSIS:  Return size of <aln>
 */
size_t ALIGNMENT_GetSize(ALIGNMENT* aln);

/*! FUNCTION:  ALIGNMENT_GetNumRegions()
 *  SYNOPSIS:  Return number of discontiguous (broken by jump state) traceback
 * alignments in <aln>. Use must have run _FindAligns() first.
 */
size_t ALIGNMENT_GetNumRegions(ALIGNMENT* aln);

/*! FUNCTION:  ALIGNMENT_GetTrace()
 *  SYNOPSIS:  Get <i>th trace in <aln>
 */
TRACE
ALIGNMENT_GetTrace(ALIGNMENT* aln, int i);

/*! FUNCTION:  ALIGNMENT_Resize()
 *  SYNOPSIS:  Resize <aln>'s trace array to <size>.
 */
void ALIGNMENT_Resize(ALIGNMENT* aln, size_t size);

/*! FUNCTION:  ALIGNMENT_GrowTo()
 *  SYNOPSIS:  Resize <aln>'s trace array to <size> only if array is smaller
 * than current size.
 */
void ALIGNMENT_GrowTo(ALIGNMENT* aln, size_t size);

/*! FUNCTION:  ALIGNMENT_Resize()
 *  SYNOPSIS:  Compare alignments <a> and <b>.
 *  RETURN:    Zero if equal.
 */
int ALIGNMENT_Compare(ALIGNMENT* a, ALIGNMENT* b);

/*! FUNCTION:  ALIGNMENT_AddTrace()
 *  SYNOPSIS:  Push trace onto end of alignment.
 */
void ALIGNMENT_AddTrace(ALIGNMENT* aln, TRACE tr);

/*! FUNCTION:  ALIGNMENT_AppendTrace()
 *  SYNOPSIS:  Append trace <st_cur, q_0, t_0> to <aln>.
 */
int ALIGNMENT_AppendTrace(ALIGNMENT* aln,   /* Traceback Alignment */
                          const int st_cur, /* HMM state */
                          const int q_0,    /* index in query/sequence */
                          const int t_0);   /* index in target/model */

/*! FUNCTION:  ALIGNMENT_AppendScore()
 *  SYNOPSIS:  Append trace <st_cur, q_0, t_0> to <aln>.
 */
int ALIGNMENT_AppendScore(ALIGNMENT* aln,     /* Traceback Alignment */
                          const float score); /* Score */

/*! FUNCTION:  ALIGNMENT_FindRegions()
 *  SYNOPSIS:  Scan full model alignment traceback for all alignement regions.
 *             These regions are those running through the core model, from a
 * BEGIN to an END state. Stores each <beg> and <end> positions in the full
 * alignment are pushed to alignment vectors <tr_beg> amd <tr_end>. Returns the
 * number of alignment regions.
 */
int ALIGNMENT_FindRegions(ALIGNMENT* aln);

/*! FUNCTION:  ALIGNMENT_ScoreRegions()
 *  SYNOPSIS:  Scan <aln> for all alignments running talo through the core model
 * (running from BEGIN to END state). Caller must run _AppendScore() for each
 * _AppendTrace().
 */
float ALIGNMENT_ScoreRegions(ALIGNMENT* aln);

/*! FUNCTION:  ALIGNMENT_AddRegion()
 *  SYNOPSIS:  Adds a distinct, discrete alignment region to list with <beg> and
 * <end> points and <score> for region.
 */
void ALIGNMENT_AddRegion(ALIGNMENT* aln, int beg, int end, float score);

/*! FUNCTION:  ALIGNMENT_SetRegion()
 *  SYNOPSIS:  Sets <beg> and <end> endpoint indexes of the <aln> alignment.
 */
void ALIGNMENT_SetRegion(ALIGNMENT* aln, int aln_idx);

/*! FUNCTION:  ALIGNMENT_Reverse()
 *  SYNOPSIS:  Reverse order of <aln>.
 *             For use with alignments built from back-to-front via backtrace.
 *             This will correct that to normal ordering.
 */
void ALIGNMENT_Reverse(ALIGNMENT* aln);

/*! FUNCTION:  ALIGNMENT_Build_HMMER_Style()
 *  SYNOPSIS:  Generate <aln> strings, HMMER-style.
 *             Stores in <target_aln>, <center_aln>, <query_aln>.
 *             Expects <aln> has already been constructed.
 */
void ALIGNMENT_Build_HMMER_Style(ALIGNMENT* aln,
                                 SEQUENCE* query,
                                 HMM_PROFILE* target);

/*! FUNCTION:  ALIGNMENT_Build_MMSEQS_Style()
 *  SYNOPSIS:  Generate <aln> strings, MMSEQS-style.
 *             Stores in <cigar_aln>.
 *             Expects <aln> has already been constructed.
 */
void ALIGNMENT_Build_MMSEQS_Style(ALIGNMENT* aln,
                                  SEQUENCE* query,
                                  HMM_PROFILE* target);

/*! FUNCTION:  ALIGNMENT_Dump()
 *  SYNOPSIS:  Outputs <aln> to open file pointer <fp>.
 */
void ALIGNMENT_Dump(ALIGNMENT* aln, FILE* fp);

/*! FUNCTION:  ALIGNMENT_Save()
 *  SYNOPSIS:  Save <aln> to file at location <filename>.
 *             Handles opening and closing of file.
 */
void ALIGNMENT_Save(ALIGNMENT* aln, char* filename);

#endif /* _ALIGNMENT_H */
