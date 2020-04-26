/*******************************************************************************
 *  FILE:      index_parser.h
 *  PURPOSE:   Build an INDEX object from: .idx, .hmm, or .fasta file
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _INDEX_PARSER_H
#define _INDEX_PARSER_H

/* index .hmm file by hmm name and offset */
F_INDEX* F_INDEX_Hmm_Build( const char* _filename_ );

/* index .hmm file by hmm name and offset */
F_INDEX* F_INDEX_Fasta_Build( const char* _filename_ );

/* load an index file */
F_INDEX* F_INDEX_Load( const char* _filename );

/* update f_index using mmseqs lookup table */
void F_INDEX_Lookup_Update( F_INDEX*   f_index, 
                            char*      _lookup_filepath_ );

/* load index file */
F_INDEX* F_INDEX_Plus_Load( const char*   _filename_ );

#endif /* _INDEX_PARSER_H */