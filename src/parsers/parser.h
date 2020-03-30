/*******************************************************************************
 *  FILE:      parsers.h
 *  PURPOSE:   Parsers for .hmm, .fasta, and .index files
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _PARSER_H
#define _PARSER_H

/* === FUNCTIONS === */

/* index .hmm file by hmm name and offset */
F_INDEX* F_INDEX_Hmm_Build(const char* _filename_);

/* index .fasta file by fasta name and offset */
F_INDEX* F_INDEX_Fasta_Build(const char* _filename_);

/* load an index file */
F_INDEX* F_INDEX_Load(const char* _filename);


SEQUENCE* SEQUENCE_Fasta_Parse(char*      _filename_, 
                               long int   offset);

#endif /* _PARSER_H */