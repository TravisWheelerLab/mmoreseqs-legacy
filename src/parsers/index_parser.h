/*******************************************************************************
 *  - FILE:      index_parser.h
 *  - DESC:    Build an INDEX object from: .idx, .hmm, or .fasta file
 *******************************************************************************/

#ifndef _INDEX_PARSER_H
#define _INDEX_PARSER_H

/*! FUNCTION:  F_INDEX_Hmm_Build()
 *  SYNOPSIS:  Build F_INDEX object from .hmm file with list of hmm name and
 * offset into file. METHOD:    [1]   scans .hmm file line-by-line [2]   if line
 * is the header of a hmm model (line starts with "HMMER") [2a]  name of hmm and
 * it offset location into a file added to index RETURN:    F_INDEX object
 * containing index of file
 */
F_INDEX* F_INDEX_Hmm_Build(F_INDEX* f_index, const char* filename);

/*! FUNCTION:  F_INDEX_Fasta_Build()
 *  SYNOPSIS:  Build F_INDEX object from .fasta file by fasta name and offset
 * into file. METHOD:    [1]   scans .fasta file line-by-line [2]   if line is
 * the header of a fasta sequence: [2a]  name of fasta and it offset location
 * into a file added to index RETURN:    F_INDEX object containing index of file
 */
F_INDEX* F_INDEX_Fasta_Build(F_INDEX* f_index, const char* filename);

/*! FUNCTION:  F_INDEX_Load()
 *  SYNOPSIS:  Build F_INDEX object from .idx file.
 *  METHOD:    [1]   scans .fasta file line-by-line
 *             [2]   line: {result_id} {offset} {}
 *  RETURN:    F_INDEX object containing index of file
 */
F_INDEX* F_INDEX_Load(F_INDEX* f_index, const char* filename);

/*! FUNCTION:  F_INDEX_Load()
 *  SYNOPSIS:  Update f_index using mmseqs lookup table
 */
void F_INDEX_Lookup_Update(F_INDEX* f_index, char* _lookup_filein_);

/*! FUNCTION:  F_INDEX_Load()
 *  SYNOPSIS:  Load F_INDEX object from .idx file.
 */
F_INDEX* F_INDEX_Plus_Load(F_INDEX* f_index, const char* filename);

#endif /* _INDEX_PARSER_H */
