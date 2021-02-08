/*******************************************************************************
 *  FILE:      m8_parser.h
 *  PURPOSE:   Parses .m8 file into RESULTS object
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _M8_PARSER_H
#define _M8_PARSER_H

/*! FUNCTION:  RESULTS_M8_Parse()
 *  SYNOPSIS:  Parse .m8 results file at <_filename_> and 
 * 				stores data in M8_RESULTS object <results>.
 * 				Only imports over the range (<start_idx>, <end_idx>).
 */
void 
RESULTS_M8_Parse( 	M8_RESULTS* 	results,
							char* 			_filename_,
							int 				start_idx,
							int 				end_idx );

#endif /* _M8_PARSER_H */