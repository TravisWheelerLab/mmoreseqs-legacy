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
RESULTS* RESULTS_M8_Parse( char* _filename_ );

#endif /* _M8_PARSER_H */