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
void SEQUENCE_Destroy(SEQUENCE *seq);

/* Set Sequence String to SEQUENCE and update length */
void SEQUENCE_Set_Seq(SEQUENCE* seq,
                      char*     seq_text);
/* Append Sequence String onto current SEQUENCE and update length */
void SEQUENCE_Append_Seq(SEQUENCE* seq,
                         char*     seq_text);

/* Set Textfield to SEQUENCE field */
void SEQUENCE_Set_Textfield(char** seq_field,
                            char*  text);

/* Convert SEQUENCE to HMM_PROFILE */
HMM_PROFILE* SEQUENCE_to_HMM_PROFILE(SEQUENCE* seq);

/* Output SEQUENCE out to FILE POINTER */
void SEQUENCE_Dump(SEQUENCE *seq,
                   FILE *fp);

#endif /* _SEQUENCE_H */
