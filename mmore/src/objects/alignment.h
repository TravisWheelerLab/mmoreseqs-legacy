/*******************************************************************************
 *  FILE:      alignment.c
 *  PURPOSE:   ALIGNMENT Object.
 *             Used for storing Viterbi or Forward/Backward Alignments.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _ALIGNMENT_H
#define _ALIGNMENT_H

/*! FUNCTION:  ALIGNMENT_Create()
 *  SYNOPSIS:  Constructs <aln>, allocates memory and returns pointer.
 */
ALIGNMENT* 
ALIGNMENT_Create();

/*! FUNCTION:  ALIGNMENT_Destroy()
 *  SYNOPSIS:  Destroys <aln>, frees memory and return NULL pointer.
 */
ALIGNMENT* 
ALIGNMENT_Destroy( ALIGNMENT*    aln );

/*! FUNCTION:  ALIGNMENT_Reuse()
 *  SYNOPSIS:  Wipes <aln>'s old data for reuse, sets dimensions <Q x T>.
 */
void 
ALIGNMENT_Reuse(     ALIGNMENT*  aln,
                     int         Q,
                     int         T );

/*! FUNCTION:  ALIGNMENT_Pushback()
 *  SYNOPSIS:  Appends <tr> to end of <aln>, resizes memory if necessary.
 */
void 
ALIGNMENT_Pushback(     ALIGNMENT*  aln,
                        TRACE*      tr );

/*! FUNCTION:  ALIGNMENT_Resize()
 *  SYNOPSIS:  Resize <aln>'s trace array to <size>.
 */
void 
ALIGNMENT_Resize(    ALIGNMENT*     aln,
                     size_t         size );

/*! FUNCTION:  ALIGNMENT_GrowTo()
 *  SYNOPSIS:  Resize <aln>'s trace array to <size> only if array is smaller than current size.
 */
void 
ALIGNMENT_GrowTo(    ALIGNMENT*   aln,
                     size_t       size );

/*! FUNCTION:  ALIGNMENT_Resize()
 *  SYNOPSIS:  Compare alignments <a> and <b>.
 *  RETURN:    Zero if equal.
 */
int ALIGNMENT_Compare(  ALIGNMENT*     a,
                        ALIGNMENT*     b );

/*! FUNCTION:  ALIGNMENT_Append()
 *  SYNOPSIS:  Append trace <st_cur, q_0, t_0> to <aln>.
 */
int 
ALIGNMENT_Append(    ALIGNMENT*   aln,       /* Traceback Alignment */
                     const int    st_cur,    /* HMM state */
                     const int    q_0,       /* index in query/sequence */
                     const int    t_0 );     /* index in target/model */

/*! FUNCTION:  ALIGNMENT_Find_Length()
 *  SYNOPSIS:  Scan <aln> for beginning, end, and length of alignment. 
 *             Stores <beg> and <end>.
 */
void 
ALIGNMENT_Find_Length(  ALIGNMENT*  aln );

/*! FUNCTION:  ALIGNMENT_SetEndpoints()
 *  SYNOPSIS:  Sets <beg> and <end> endpoint indexes of the <aln> alignment.
 */
void 
ALIGNMENT_SetEndpoints(   ALIGNMENT*  aln,
                           int         beg,
                           int         end );

/*! FUNCTION:  ALIGNMENT_Reverse()
 *  SYNOPSIS:  Reverse order of <aln>.
 *             For use with alignments built from back-to-front via backtrace. 
 *             This will correct that to normal ordering.
 */
void 
ALIGNMENT_Reverse( ALIGNMENT*    aln );

/*! FUNCTION:  ALIGNMENT_Build_HMMER_Style()
 *  SYNOPSIS:  Generate <aln> strings, HMMER-style. 
 *             Stores in <target_aln>, <center_aln>, <query_aln>.
 *             Expects <aln> has already been constructed.
 */
void 
ALIGNMENT_Build_HMMER_Style(     ALIGNMENT*     aln,
                                 SEQUENCE*      query,
                                 HMM_PROFILE*   target );

/*! FUNCTION:  ALIGNMENT_Build_MMSEQS_Style()
 *  SYNOPSIS:  Generate <aln> strings, MMSEQS-style. 
 *             Stores in <cigar_aln>.
 *             Expects <aln> has already been constructed.
 */
void 
ALIGNMENT_Build_MMSEQS_Style(    ALIGNMENT*     aln,
                                 SEQUENCE*      query,
                                 HMM_PROFILE*   target );

/*! FUNCTION:  ALIGNMENT_Dump()
 *  SYNOPSIS:  Outputs <aln> to open file pointer <fp>.
 */
void 
ALIGNMENT_Dump(   ALIGNMENT*  aln,
                  FILE*       fp );

/*! FUNCTION:  ALIGNMENT_Save()
 *  SYNOPSIS:  Save <aln> to file at location <_filename_>. 
 *             Handles opening and closing of file.
 */
void 
ALIGNMENT_Save(   ALIGNMENT*  aln,
                  char*       _filename_ );

#endif /* _ALIGNMENT_H */