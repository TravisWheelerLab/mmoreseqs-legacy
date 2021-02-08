/*******************************************************************************
 *  FILE:      sequence.h
 *  PURPOSE:   SEQUENCE object
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _SEQUENCE_H
#define _SEQUENCE_H

/*! FUNCTION:  SEQUENCE_Create()
 *  SYNOPSIS:  Create <seq>, allocate memory and return pointer.
 */
SEQUENCE* 
SEQUENCE_Create();

/*! FUNCTION:  SEQUENCE_Destroy()
 *  SYNOPSIS:  Destroy <seq>, free memory and return NULL pointer.
 */
SEQUENCE* 
SEQUENCE_Destroy( SEQUENCE*  seq );

/*! FUNCTION:  SEQUENCE_Reuse()
 *  SYNOPSIS:  Reuse <seq> by reinitializing all fields except <seq> field.
 */
STATUS_FLAG
SEQUENCE_Reuse( SEQUENCE* seq );

/** FUNCTION:  SEQUENCE_GetSeq()
 *  SYNOPSIS:  Get Sequence String <seq>.
 */
STR 
SEQUENCE_GetSeq(  SEQUENCE*   seq );

/* Set Sequence String to SEQUENCE and update length */
void 
SEQUENCE_SetSeq( SEQUENCE*  seq,
                  char*      seq_text );

/* Append Sequence String onto current SEQUENCE and update length */
void 
SEQUENCE_Append_Seq( SEQUENCE*  seq,
                     char*      seq_text);

/* Reallocate space for SEQUENCE */
void 
SEQUENCE_Resize_Seq( SEQUENCE*    seq,
                     int          size );

/* Set Textfield to SEQUENCE field */
void 
SEQUENCE_SetTextfield( STR*     seq_field,
                        STR      text);

/* Set sequence to a subsequence */
void 
SEQUENCE_SetDomain( SEQUENCE*   seq, 
                     RANGE       Q_range );

/* Unset sequence to domain */
void 
SEQUENCE_UnsetDomain( SEQUENCE*  seq );

/* Output SEQUENCE out to FILE POINTER */
void 
SEQUENCE_Dump( SEQUENCE*  seq,
               FILE*      fp);

#endif /* _SEQUENCE_H */
