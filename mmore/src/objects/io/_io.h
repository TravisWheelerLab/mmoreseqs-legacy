/*******************************************************************************
 *  FILE:      _io.h
 *  PURPOSE:   Basic input/output objects (readers, writers).
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _IO_H
#define _IO_H

/* filer for file management (open/close/etc) */
#include "filer.h"
/* buffer for storing file data */
#include "buffer.h"
#include "buffer_read.h"
#include "buffer_write.h"

/* reader and writer (depend on buffer) */
#include "reader.h"
#include "writer.h"

#endif /* _IO_H */
