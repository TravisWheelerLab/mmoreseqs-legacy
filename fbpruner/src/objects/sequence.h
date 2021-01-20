/*******************************************************************************
 *  FILE:      sequence.h
 *  PURPOSE:   SEQUENCE object
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _SEQUENCE_H
#define _SEQUENCE_H

/** FUNCTION:  SEQUENCE_Create()
 *  SYNOPSIS:  
 *
 *  RETURN:    Return pointer to new SEQUENCE object.
 */
SEQUENCE* 
SEQUENCE_Create();

/** FUNCTION:  SEQUENCE_Destroy()
 *  SYNOPSIS:  
 *
 *  RETURN:    Return NULL pointer.
 */
SEQUENCE* 
SEQUENCE_Destroy( SEQUENCE*  seq );

/* Reuse sequence by reinitializing all fields except seq field */
void 
SEQUENCE_Reuse( SEQUENCE* seq );

/* Set Sequence String to SEQUENCE and update length */
void 
SEQUENCE_Set_Seq( SEQUENCE*  seq,
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
SEQUENCE_Set_Textfield( char**  seq_field,
                        char*   text);

/* Set sequence to a subsequence */
void 
SEQUENCE_Set_Domain( SEQUENCE*   seq, 
                     RANGE       Q_range );

void 
SEQUENCE_Unset_Domain( SEQUENCE*  seq );

/* Output SEQUENCE out to FILE POINTER */
void 
SEQUENCE_Dump( SEQUENCE*  seq,
               FILE*      fp);

#endif /* _SEQUENCE_H */
