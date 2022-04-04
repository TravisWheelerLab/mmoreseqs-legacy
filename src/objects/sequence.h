/*******************************************************************************
 *  - FILE:      sequence.h
 *  - DESC:    SEQUENCE object
 *******************************************************************************/

#ifndef _SEQUENCE_H
#define _SEQUENCE_H

/*! FUNCTION:  SEQUENCE_Create()
 *  SYNOPSIS:  Create <seq>, allocate memory and return pointer.
 */
SEQUENCE* SEQUENCE_Create();

/*! FUNCTION:  SEQUENCE_Destroy()
 *  SYNOPSIS:  Destroy <seq>, free memory and return NULL pointer.
 */
SEQUENCE* SEQUENCE_Destroy(SEQUENCE* seq);

/*! FUNCTION:  SEQUENCE_Reuse()
 *  SYNOPSIS:  Reuse <seq> by reinitializing all fields except <seq> field.
 */
STATUS_FLAG
SEQUENCE_Reuse(SEQUENCE* seq);

/*! FUNCTION:  SEQUENCE_GetSize()
 *  SYNOPSIS:  Gets the length of the sequence <seq>.
 */
size_t SEQUENCE_GetSize(SEQUENCE* seq);

/** FUNCTION:  SEQUENCE_GetSeq()
 *  SYNOPSIS:  Get Sequence String <seq>.
 */
STR SEQUENCE_GetSeq(SEQUENCE* seq);

/* Set Sequence String to SEQUENCE and update length */
void SEQUENCE_SetSeq(SEQUENCE* seq, char* seq_text);

/* Append Sequence String onto current SEQUENCE and update length */
void SEQUENCE_AppendSeq(SEQUENCE* seq, char* seq_text);

/* Reallocate space for SEQUENCE */
void SEQUENCE_Resize(SEQUENCE* seq, int size);

/** FUNCTION:  SEQUENCE_GetCharAt()
 *  SYNOPSIS:  Get character <c> from sequence <seq> at <i>th position.
 */
char SEQUENCE_GetCharAt(const SEQUENCE* seq, const int i);

/** FUNCTION:  SEQUENCE_GetDigitAt()
 *  SYNOPSIS:  Get int <val> from digitized sequence at <i>th position.
 *             Caller must call _Digitize() before this.
 */
int SEQUENCE_GetDigitAt(const SEQUENCE* seq, const int i);

/** FUNCTION:  SEQUENCE_Resize()
 *  SYNOPSIS:  Set Textfield <test> to SEQUENCE field <seq_field> (overwrites).
 *
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors.
 */
void SEQUENCE_SetTextfield(STR* seq_field, STR text);

/** FUNCTION:  SEQUENCE_Digitize()
 *  SYNOPSIS:  Digitize text <seq> to create digital sequence <dseq>.
 *
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors.
 */
void SEQUENCE_Digitize(SEQUENCE* seq);

/** FUNCTION:  SEQUENCE_SetDomain()
 *  SYNOPSIS:  Set SEQUENCE <seq> to cover a subsequence <seq> to sequence
 * <full_seq>. Subsequence covers domain range <q_beg, q_end>.
 *
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors.
 */
void SEQUENCE_SetDomain(SEQUENCE* seq, RANGE Q_range);

/** FUNCTION:  SEQUENCE_UnsetSubseq()
 *  SYNOPSIS:  Set SEQUENCE <seq> to cover a subsequence <seq> to sequence
 * <full_seq>. Subsequence covers range <q_beg, q_end>.
 *
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors.
 */
void SEQUENCE_UnsetDomain(SEQUENCE* seq);

/** FUNCTION:  SEQUENCE_Dump()
 *  SYNOPSIS:  Set SEQUENCE <seq> to file pointer <fp>.
 *
 *  RETURN:    Returns <STATUS_SUCCESS> if no errors.
 */
void SEQUENCE_Dump(SEQUENCE* seq, FILE* fp);

#endif /* _SEQUENCE_H */
