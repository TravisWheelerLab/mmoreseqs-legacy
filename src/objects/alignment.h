/*******************************************************************************
 *  FILE:    alignment.c
 *  PURPOSE: ALIGNMENT Object.
 *
 *  AUTHOR:  Dave Rich
 *  BUG:     Lots.
 *******************************************************************************/

#ifndef _ALIGNMENT_H
#define _ALIGNMENT_H

/* === OBJECTS === */
// typedef struct {
//    int         i;            /* index in query */
//    int         j;            /* index in target */
//    int         st;           /* state at index */
// } TRACE;

// typedef struct {
//    int         N;            /* current length */
//    int         Nalloc;       /* allocated length */
//    int         beg;          /* position in trace for first MID state */
//    int         end;          /* position in trace for last MID state */
//    TRACE       *traces;      /* list of all (state,i,j) TRACES in ALIGNMENT */
// } ALIGNMENT;

/* === FUNCTIONS === */
/* constructor */
ALIGNMENT* ALIGNMENT_Create();
/* destructor */
void ALIGNMENT_Destroy(ALIGNMENT *aln);

/* push trace onto end of alignment */
void ALIGNMENT_Pushback(ALIGNMENT* aln,
                        TRACE*     tr);
/* resize TRACE array in ALIGNMENT */
void ALIGNMENT_Resize(ALIGNMENT* aln,
                      float      growth_factor);

/* outputs ALIGNMENT to FILE pointer */
void ALIGNMENT_Dump(ALIGNMENT* aln,
                    FILE*      fp);
/* saves ALIGNMENT to FILE with filename */
void ALIGNMENT_Save(ALIGNMENT* aln,
                    char*      _filename_);

#endif /* _ALIGNMENT_H */