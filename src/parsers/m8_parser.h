/*******************************************************************************
 *  - FILE:   m8_parser.h
 *  - DESC:    Parses .m8 file into RESULTS object
 *******************************************************************************/

#ifndef _M8_PARSER_H
#define _M8_PARSER_H

/*! FUNCTION:  RESULTS_M8_Parse()
 *  SYNOPSIS:  Parse .m8 results file at <filename> and
 * 				stores data in M8_RESULTS object <results>.
 * 				Only imports over the range (<start_idx>, <end_idx>).
 */
void RESULTS_M8_Parse(M8_RESULTS* results,
                      char* filename,
                      int start_idx,
                      int end_idx);

#endif /* _M8_PARSER_H */
