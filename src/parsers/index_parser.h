/*******************************************************************************
 *  FILE:      index_parser.h
 *  PURPOSE:   Build an INDEX object from: .idx, .hmm, or .fasta file
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _INDEX_PARSER_H
#define _INDEX_PARSER_H

/* === OBJECTS === */

// typedef struct {
//    int      N;
//    int      Nalloc;
//    long*    offsets;
//    char**   names;
// } INDEX_HMM;

/* === FUNCTIONS === */

/* index .hmm file by hmm name and offset */
F_INDEX* F_INDEX_Hmm_Build( const char* _filename_ );

/* index .hmm file by hmm name and offset */
F_INDEX* F_INDEX_Fasta_Build( const char* _filename_ );

/* load an index file */
F_INDEX* F_INDEX_Load( const char* _filename );

#endif /* _INDEX_PARSER_H */