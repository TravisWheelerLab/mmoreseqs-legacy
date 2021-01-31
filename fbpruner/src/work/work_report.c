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
   ARGS*    args  = worker->args;

   /* TODO: need to handle results in bulk. currently report one-at-a-time. */
   /* open file pointers */
   if ( args->is_redirect_stdout ) {
      worker->output_fp = ERROR_fopen( args->output_filepath, "w" );
   } else {
      worker->output_fp = stdout;
   }
   if ( args->is_tblout ) {
      worker->tblout_fp = ERROR_fopen( args->tblout_filepath, "w" );
   }
   if ( args->is_m8out ) {
      worker->m8out_fp = ERROR_fopen( args->m8out_filepath, "w" );
   }
   if ( args->is_myout ) {
      worker->myout_fp = ERROR_fopen( args->myout_filepath, "w" );
   }
   if ( args->is_mydomout ) {
      worker->mydomout_fp  = ERROR_fopen( args->mydomout_filepath, "w" );
   }
   if ( args->is_mytimeout ) {
      worker->mytimeout_fp = ERROR_fopen( args->mytimeout_filepath, "w" );
   }
   if ( args->is_mythreshout ) {
      worker->mythreshout_fp = ERROR_fopen( args->mythreshout_filepath, "w" );
   }
}

/*! FUNCTION:  	WORK_close()
 *  SYNOPSIS:  	Close all open files in <worker>.
 */
void 
WORK_close( WORKER* worker )
{
   ARGS*    args  =  worker->args;  

   if ( args->is_redirect_stdout ) {
      if ( worker->output_fp != stdout ) {
         worker->output_fp    = ERROR_fclose( worker->output_fp );
      } 
   }
   if ( args->is_tblout ) {
      worker->tblout_fp    = ERROR_fclose( worker->tblout_fp );
   }
   if ( args->is_m8out ) {
      worker->m8out_fp     = ERROR_fclose( worker->m8out_fp );
   }
   if ( args->is_myout ) {
      worker->myout_fp    = ERROR_fclose( worker->myout_fp );
   }
   if ( args->is_mydomout ) {
      worker->mydomout_fp = ERROR_fclose( worker->mydomout_fp );
   }
   if ( args->is_mytimeout ) {
      worker->mytimeout_fp = ERROR_fclose( worker->mytimeout_fp );
   }
   if ( args->is_mythreshout ) {
      worker->mythreshout_fp = ERROR_fclose( worker->mythreshout_fp );
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
   REPORT_stdout_header( worker, worker->output_fp );

   if ( args->is_tblout ) {
      REPORT_domtblout_header( worker, worker->tblout_fp );
   }
   if ( args->is_m8out ) {
      REPORT_m8out_header( worker, worker->m8out_fp );
   }
   if ( args->is_myout ) {
      REPORT_myout_header( worker, worker->myout_fp );
   }
   if ( args->is_mytimeout ) {
      REPORT_mytimeout_header( worker, worker->mytimeout_fp );
   }
   if ( args->is_mythreshout ) {
      REPORT_mythreshout_header( worker, worker->mythreshout_fp );
   }
   if ( args->is_mydomout && args->is_run_domains ) {
      REPORT_mydomout_header( worker, worker->mydomout_fp );
   }
}

/*! FUNCTION:  	WORK_report_result_current()
 *  SYNOPSIS:  	Write current result entry to all open files in <worker>.
 */
void 
WORK_report_result_current( WORKER*  worker )
{
   ARGS* args = worker->args;

   /* print result to all open pointers */
   REPORT_stdout_entry( worker, worker->result, worker->output_fp );

   if ( args->is_tblout ) {
      REPORT_domtblout_entry( worker, worker->result, worker->tblout_fp );
   }
   if ( args->is_m8out ) {
      REPORT_m8out_entry( worker, worker->result, worker->m8out_fp );
   }
   if ( args->is_myout ) { 
      REPORT_myout_entry( worker, worker->result, worker->myout_fp );
   }
   if ( args->is_mytimeout ) { 
      REPORT_mytimeout_entry( worker, worker->result, worker->mytimeout_fp );
   }
   if ( args->is_mythreshout ) { 
      REPORT_mythreshout_entry( worker, worker->result, worker->mythreshout_fp );
   }
   if ( args->is_mydomout && args->is_run_domains ) { 
      REPORT_mydomout_entry( worker, worker->result, worker->mydomout_fp );
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
   REPORT_stdout_footer( worker, worker->output_fp );

   if ( args->is_tblout ) {
      REPORT_domtblout_footer( worker, worker->tblout_fp );
   }
   if ( args->is_m8out ) {
      REPORT_m8out_footer( worker, worker->m8out_fp );
   }
   if ( args->is_myout ) { 
      REPORT_myout_footer( worker, worker->myout_fp );
   }
   if ( args->is_mytimeout ) { 
      REPORT_mytimeout_footer( worker, worker->mytimeout_fp );
   }
   if ( args->is_mythreshout ) { 
      REPORT_mythreshout_footer( worker, worker->mythreshout_fp );
   }
   if ( args->is_mydomout && args->is_run_domains ) { 
      REPORT_myout_footer( worker, worker->mydomout_fp );
   }
}