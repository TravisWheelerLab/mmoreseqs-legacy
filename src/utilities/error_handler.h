/*******************************************************************************
 *  FILE:      error_handler.h
 *  PURPOSE:   Functions for handling error calls
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _ERROR_HANDLER_H
#define _ERROR_HANDLER_H

/* outputs error to console and terminates program */
void ERROR_Handler( const int    error_code,
                    const char*  file,
                    const int    line,
                    const char*  func,
                    const char*  note );

void ERROR_NullPtrCheck( const void*   data,
                         int*          error_code,
                         const char*   note );

#endif /* _ERROR_HANDLER_H */