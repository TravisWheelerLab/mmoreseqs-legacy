/*******************************************************************************
 *  FILE:      work_report.c
 *  PURPOSE:   Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Subroutines for outputting reports.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

/* local imports */
#include "../objects/structs.h"
#include "../utilities/_utilities.h"
#include "../objects/_objects.h"
#include "../parsers/_parsers.h"
#include "../algs_linear/_algs_linear.h"
#include "../algs_quad/_algs_quad.h"
#include "../algs_naive/_algs_naive.h"
#include "../algs_sparse/_algs_sparse.h"
#include "../reporting/_reporting.h"

/* headers */
#include "_work.h"
#include "work_report.h"

/*! FUNCTION:  	WORK_open()
 *  SYNOPSIS:  	Open all valid files in <worker>.
 */
void 
WORK_open( WORKER* worker )
{
   ARGS* args = worker->args; 

   /* open file pointers */
   if ( args->is_redirect_stdout ) {
      FILER_Open( worker->output_file );
   }
   if ( worker->tblout_file != NULL && args->is_tblout ) {
      FILER_Open( worker->tblout_file );
   }
   if ( worker->m8out_file != NULL && args->is_m8out ) {
      FILER_Open( worker->m8out_file );
   }
   if ( worker->myout_file != NULL && args->is_myout ) {
      FILER_Open( worker->myout_file );
   }
   if ( worker->mydomout_file != NULL && args->is_mydomout ) {
      FILER_Open( worker->mydomout_file );
   }
   if ( worker->mytimeout_file != NULL && args->is_mytimeout ) {
      FILER_Open( worker->mytimeout_file );
   }
   if ( worker->mythreshout_file != NULL && args->is_mythreshout ) {
      FILER_Open( worker->mythreshout_file );
   }
}

/*! FUNCTION:  	WORK_close()
 *  SYNOPSIS:  	Close all open files in <worker>.
 */
void 
WORK_close( WORKER* worker )
{
   ARGS* args = worker->args; 

   /* open file pointers */
   if ( args->is_redirect_stdout ) {
      FILER_Close( worker->output_file );
   }
   if ( args->is_tblout ) {
      FILER_Close( worker->tblout_file );
   }
   if ( args->is_m8out ) {
      FILER_Close( worker->m8out_file );
   }
   if ( args->is_myout ) {
      FILER_Close( worker->myout_file );
   }
   if ( args->is_mydomout ) {
      FILER_Close( worker->mydomout_file );
   }
   if ( args->is_mytimeout ) {
      FILER_Close( worker->mytimeout_file );
   }
   if ( args->is_mythreshout ) {
      FILER_Close( worker->mythreshout_file );
   }
}

/*! FUNCTION:  	WORK_report_header()
 *  SYNOPSIS:  	Write all report headers to all open files in <worker>. 
 */
void 
WORK_report_header( WORKER*    worker )
{
   ARGS* args = worker->args;

   /* print footers to all open pointers */
   REPORT_stdout_header( worker, worker->output_file->fp );

   if ( args->is_tblout ) {
      REPORT_domtblout_header( worker, worker->tblout_file->fp );
   }
   if ( args->is_m8out ) {
      REPORT_m8out_header( worker, worker->m8out_file->fp );
   }
   if ( args->is_myout ) {
      REPORT_myout_header( worker, worker->myout_file->fp );
   }
   if ( args->is_mytimeout ) {
      REPORT_mytimeout_header( worker, worker->mytimeout_file->fp );
   }
   if ( args->is_mythreshout ) {
      REPORT_mythreshout_header( worker, worker->mythreshout_file->fp );
   }
   if ( args->is_mydomout && args->is_run_domains ) {
      REPORT_mydomout_header( worker, worker->mydomout_file->fp );
   }
}

/*! FUNCTION:  	WORK_report_result_current()
 *  SYNOPSIS:  	Write current result entry to all open files in <worker>.
 */
void 
WORK_report_result_current( WORKER*  worker )
{
   ARGS*    args     = worker->args;
   RESULT*  result   = worker->result;

   /* only add entry to these reports if search passed reporting threshold */
   if ( result->is_passed_fwdback == true ) 
   {
      REPORT_stdout_entry( worker, worker->result, worker->output_file->fp );

      if ( args->is_tblout ) {
         REPORT_domtblout_entry( worker, worker->result, worker->tblout_file->fp );
      }
      if ( args->is_m8out ) {
         REPORT_m8out_entry( worker, worker->result, worker->m8out_file->fp );
      }
      if ( args->is_myout ) { 
         REPORT_myout_entry( worker, worker->result, worker->myout_file->fp );
      }
      if ( args->is_mydomout && args->is_run_domains ) { 
         REPORT_mydomout_entry( worker, worker->result, worker->mydomout_file->fp );
      }
   }

   if ( args->is_mytimeout ) { 
         REPORT_mytimeout_entry( worker, worker->result, worker->mytimeout_file->fp );
      }
   if ( args->is_mythreshout ) { 
      REPORT_mythreshout_entry( worker, worker->result, worker->mythreshout_file->fp );
   }
}

/*! FUNCTION:  	WORK_report_result_all()
 *  SYNOPSIS:  	Write all result entries in queue to all open files in <worker>.
 */
void 
WORK_report_result_all( WORKER*   worker )
{
   ARGS* args = worker->args;
}

/*! FUNCTION:  	WORK_report_footer()
 *  SYNOPSIS:  	Write all report footers to all open files in <worker>. 
 */
void 
WORK_report_footer( WORKER*    worker )
{
   ARGS* args = worker->args;

   /* print footers to all open pointers */
   REPORT_stdout_footer( worker, worker->output_file->fp );

   if ( args->is_tblout ) {
      REPORT_domtblout_footer( worker, worker->tblout_file->fp );
   }
   if ( args->is_m8out ) {
      REPORT_m8out_footer( worker, worker->m8out_file->fp );
   }
   if ( args->is_myout ) { 
      REPORT_myout_footer( worker, worker->myout_file->fp );
   }
   if ( args->is_mytimeout ) { 
      REPORT_mytimeout_footer( worker, worker->mytimeout_file->fp );
   }
   if ( args->is_mythreshout ) { 
      REPORT_mythreshout_footer( worker, worker->mythreshout_file->fp );
   }
   if ( args->is_mydomout && args->is_run_domains ) { 
      REPORT_myout_footer( worker, worker->mydomout_file->fp );
   }
}