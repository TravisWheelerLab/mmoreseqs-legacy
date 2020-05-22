/*******************************************************************************
 *  FILE:      seq_parser.h
 *  PURPOSE:   Parses .fasta files into SEQUENCE object
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _SEQ_PARSER_H
#define _SEQ_PARSER_H

/* === INCLUDES === */
// #include "objects/structs.c"
// #include "objects/sequence.c"
// #include "objects/hmm_profile.c"

/* === FUNCTIONS === */
/* parse .fasta file and build SEQUENCE object */
void SEQUENCE_Fasta_Parse( SEQUENCE*   seq,
                           char*       _filename_, 
                           long        offset );

#endif /* _SEQ_PARSER_H */
