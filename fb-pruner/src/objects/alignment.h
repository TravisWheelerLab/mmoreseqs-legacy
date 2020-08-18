/*******************************************************************************
 *  FILE:      alignment.c
 *  PURPOSE:   ALIGNMENT Object.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:
 *******************************************************************************/

#ifndef _ALIGNMENT_H
#define _ALIGNMENT_H

/* constructor */
ALIGNMENT* ALIGNMENT_Create();

/* destructor */
void* ALIGNMENT_Destroy( ALIGNMENT*    aln );

/* reuse alignment by clearing traceback and setting new dimensions */
void ALIGNMENT_Reuse(   ALIGNMENT*  aln,
                        int         Q,
                        int         T );

/* push trace onto end of alignment */
void ALIGNMENT_Pushback(   ALIGNMENT* aln,
                           TRACE*     tr );

/* resize TRACE array in ALIGNMENT */
void ALIGNMENT_Resize(  ALIGNMENT*     aln,
                        int            size );

/* resize TRACE array in ALIGNMENT */
void ALIGNMENT_GrowTo( ALIGNMENT*   aln,
                       int          size );

/* compare two alignments */
int ALIGNMENT_Compare(  ALIGNMENT*     a,
                        ALIGNMENT*     b );

/*
 *  FUNCTION:  traceback_Append()
 *  SYNOPSIS:  Append next state to Optimal Alignment.
 *  RETURN:    Return <STATUS_SUCCESS> if no errors.
 */
int ALIGNMENT_Append(   ALIGNMENT*   aln,       /* Traceback Alignment */
                        TRACE*       tr,        /* Traceback being Appended */
                        const int    st_cur,    /* HMM state */
                        const int    q_0,       /* index in query/sequence */
                        const int    t_0 );     /* index in target/model */

/* Reverse order of ALIGNMENT */
void ALIGNMENT_Reverse( ALIGNMENT*  aln );

/* create string of alignment */
void ALIGNMENT_Build_String(  ALIGNMENT*     aln,
                              SEQUENCE*      query,
                              HMM_PROFILE*   target );

/* outputs ALIGNMENT to FILE pointer */
void ALIGNMENT_Dump(    ALIGNMENT*  aln,
                        FILE*       fp );

/* saves ALIGNMENT to FILE with filename */
void ALIGNMENT_Save(    ALIGNMENT*  aln,
                        char*       _filename_ );

#endif /* _ALIGNMENT_H */