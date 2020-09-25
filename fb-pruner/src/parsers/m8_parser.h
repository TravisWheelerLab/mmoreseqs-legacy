/*******************************************************************************
 *  FILE:      m8_parser.h
 *  PURPOSE:   Parses .m8 file into RESULTS object
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _M8_PARSER_H
#define _M8_PARSER_H

/* parse .m8 results file and create RESULTS object */
void RESULTS_M8_Parse( RESULTS* 	results,
					   char* 		_filename_,
					   int 			start_idx,
					   int 			end_idx );

/* */
void RESULTS_M8_Plus_Parse( RESULTS* 	results,
					   		char* 		_filename_ );

#endif /* _M8_PARSER_H */