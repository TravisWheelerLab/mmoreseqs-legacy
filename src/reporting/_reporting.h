/*******************************************************************************
 *  - FILE:  _reporting.h
 *  - DESC:  Tools for generating and outputting reports.
 *******************************************************************************/

#ifndef _REPORTING_H
#define _REPORTING_H

/* utility functions used by all formats */
#include "report_util.h"
/* Main output */
#include "mainout.h"
/* HMMER-style output formats */
#include "domtblout.h"
#include "hmmerout.h"
/* BLAST / MMSEQS-style output formats */
#include "m8out.h"
/* Custom-style output formats */
#include "mydomout.h"
#include "myout.h"
#include "mythreshout.h"
#include "mytimeout.h"

#endif /* _REPORTING_H */
