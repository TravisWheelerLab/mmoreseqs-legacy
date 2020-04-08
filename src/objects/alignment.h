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
void ALIGNMENT_Destroy(ALIGNMENT *aln);

/* push trace onto end of alignment */
void ALIGNMENT_Pushback(ALIGNMENT* aln,
                        TRACE*     tr);
/* resize TRACE array in ALIGNMENT */
void ALIGNMENT_Resize(ALIGNMENT* aln,
                      int        size);

/* Empty ALIGNMENT Array */
void ALIGNMENT_Clear(ALIGNMENT* aln);

/* outputs ALIGNMENT to FILE pointer */
void ALIGNMENT_Dump(ALIGNMENT* aln,
                    FILE*      fp);
/* saves ALIGNMENT to FILE with filename */
void ALIGNMENT_Save(ALIGNMENT* aln,
                    char*      _filename_);

#endif /* _ALIGNMENT_H */