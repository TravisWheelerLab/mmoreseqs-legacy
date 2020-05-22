/*******************************************************************************
 *  @file sequence.c
 *  @brief SEQUENCE Object
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

#ifndef _SEQUENCE_H
#define _SEQUENCE_H

/* === OBJECTS === */
// typedef struct {
//    int    N;
//    char*  filename;
//    char*  name;
//    char*  alph;
//    char*  seq;
// } SEQUENCE;

/* === FUNCTIONS === */

/* Constructor */
SEQUENCE* SEQUENCE_Create();

/* Destructor */
void SEQUENCE_Destroy(SEQUENCE*  seq);

/* Reuse sequence by reinitializing all fields except seq field */
void SEQUENCE_Reuse(SEQUENCE* seq);

/* Set Sequence String to SEQUENCE and update length */
void SEQUENCE_Set_Seq(SEQUENCE*  seq,
                      char*      seq_text);

/* Append Sequence String onto current SEQUENCE and update length */
void SEQUENCE_Append_Seq(SEQUENCE*  seq,
                         char*      seq_text);

/* Reallocate space for SEQUENCE */
void SEQUENCE_Resize_Seq( SEQUENCE*    seq,
                          int          size );

/* Set Textfield to SEQUENCE field */
void SEQUENCE_Set_Textfield(char**  seq_field,
                            char*   text);

/* Output SEQUENCE out to FILE POINTER */
void SEQUENCE_Dump(SEQUENCE*  seq,
                   FILE*      fp);

#endif /* _SEQUENCE_H */