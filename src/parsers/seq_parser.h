/*******************************************************************************
 *  - FILE:   seq_parser.h
 *  - DESC:    Parses .fasta files into SEQUENCE object
 *******************************************************************************/

#ifndef _SEQ_PARSER_H
#define _SEQ_PARSER_H

/* parse .fasta file and build SEQUENCE object */
void SEQUENCE_Fasta_Parse(SEQUENCE* seq, char* filename, long offset);

#endif /* _SEQ_PARSER_H */
