/*******************************************************************************
 *  FILE:      _reporting.h
 *  PURPOSE:   Tools for generating and outputting reports.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _REPORTING_H
#define _REPORTING_H

/* utility functions used by all formats */
#include "report_util.h"
/* HMMER-style output formats */
#include "mainout.h"
#include "domtblout.h"
/* BLAST / MMSEQS-style output formats */
#include "m8out.h"
/* Custom-style output formats */
#include "myout.h"
#include "mydomout.h"
#include "mytimeout.h"
#include "mythreshout.h"

#endif /* _REPORTING_H */