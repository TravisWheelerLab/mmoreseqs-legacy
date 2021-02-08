/*******************************************************************************
 *  FILE:      work_etc.c
 *  PURPOSE:   Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various functions.
 *             Various uncategorized functions.
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

/* header */
#include "_work.h"
#include "work_etc.h"

/*! FUNCTION:  	WORK_set_ranges()
 *  SYNOPSIS:  	Set and Verify <worker> ranges based on <args>.
 *                If <args> query and target ranges fall out-of-bounds, constrains range to valid range.
 */
void 
WORK_set_ranges( WORKER*    worker )
{
   ARGS* args = worker->args;

   /* verify range of targets to search */
   /* if beginning is less than zero, default to entire range */
   if ( args->t_range.beg < 0 ) {
      args->t_range = (RANGE) { 0, worker->t_index->N };
   } 
   /* check for invalid ranges ( assert min > max, and min/max are less than the total in file ) */
   else if ( args->t_range.beg > args->t_range.end  || 
           args->t_range.beg > worker->t_index->N || 
           args->t_range.end > worker->t_index->N ) {
      fprintf(stderr, "ERROR: Invalid target id range (%d,%d).\n", args->t_range.beg, args->t_range.end );
      exit(EXIT_FAILURE);
   }

   /* verify range of queries to search */
   /* if beginning is less than zero, default to entire range */
   if ( args->q_range.beg < 0 ) {
      args->q_range = (RANGE) { 0, worker->q_index->N };
   } 
   /* check for invalid ranges ( assert min > max, and min/max are less than the total in file ) */
   else if ( args->q_range.beg > args->q_range.end  || 
           args->q_range.beg > worker->q_index->N || 
           args->q_range.end > worker->q_index->N ) {
      fprintf(stderr, "ERROR: Invalid target id range (%d,%d).\n", args->q_range.beg, args->q_range.end );
      exit(EXIT_FAILURE);
   }
}

/* initialize all scores to -inf */
void WORK_scores_init(  WORKER*        worker,
                        ALL_SCORES*    scores )
{
   const float val = -INF;

   /* naive algs */
   scores->naive_bound_fwd          = val;
   scores->naive_bound_bck          = val;
   /* quadratic algs */
   scores->quad_fwd                 = val; 
   scores->quad_bck                 = val;
   scores->quad_vit                 = val; 
   scores->quad_cloud_fwd           = val;
   scores->quad_cloud_bck           = val;
   scores->quad_bound_fwd           = val;
   scores->quad_bound_bck           = val;
   /* linear algs */   
   scores->lin_fwd                  = val;
   scores->lin_bck                  = val;
   scores->lin_vit                  = val; 
   scores->lin_cloud_fwd            = val;
   scores->lin_cloud_bck            = val;
   scores->lin_bound_fwd            = val;
   scores->lin_bound_bck            = val;
   /* sparse algs */
   scores->sparse_fwd               = val;
   scores->sparse_bck               = val;
   scores->sparse_vit               = val; 
   scores->sparse_bound_fwd         = val;
   scores->sparse_bound_bck         = val;
   /* threshold scores */
   scores->threshold_vit            = val; 
   scores->threshold_cloud_max      = val;
   scores->threshold_cloud_compo    = val;
   scores->threshold_bound_max      = val;
   scores->threshold_dom_max        = val;  
   scores->threshold_dom_compo      = val;   
}

/* initialize all times to zero */
void 
WORK_times_init(  WORKER*    worker,
                  TIMES*     times )
{
   const float val = 0.0f;

   /* total */
   times->program_start        = val;
   times->program_end          = val;
   times->program              = val;
   /* indexes */
   times->load_target_index    = val;
   times->load_query_index     = val;
   /* main loop */
   times->loop_start           = val;
   times->loop_end             = val;
   times->loop                 = val;
   /* load */
   times->load_target          = val;
   times->load_query           = val;
   /* naive algs */
   times->naive_cloud          = val;        
   /* quadratic algs */
   times->quad_fwd             = val;
   times->quad_bck             = val;
   times->quad_vit             = val;
   times->quad_trace           = val;
   times->quad_cloud_fwd       = val;
   times->quad_cloud_bck       = val;
   times->quad_merge           = val;
   times->quad_reorient        = val;
   times->quad_bound_fwd       = val;
   times->quad_bound_bck       = val;
   times->quad_posterior       = val;
   times->quad_optacc          = val;
   /* linear algs */   
   times->lin_fwd              = val;
   times->lin_bck              = val;
   times->lin_vit              = val;
   times->lin_trace            = val;
   times->lin_cloud_fwd        = val;
   times->lin_cloud_bck        = val;
   times->lin_merge            = val;
   times->lin_reorient         = val;
   times->lin_bound_fwd        = val;
   times->lin_bound_bck        = val;
   /* sparse algs */
   times->sp_build_mx          = val;
   times->sp_fwd               = val;
   times->sp_bck               = val;
   times->sp_vit               = val;
   times->sp_trace             = val;
   times->sp_cloud_fwd         = val;
   times->sp_cloud_bck         = val;
   times->sp_bound_fwd         = val;
   times->sp_bound_bck         = val;
   times->sp_posterior         = val;
   times->sp_optacc            = val;
   /* doms */
   times->dom_total            = val;
   times->dom_bound_fwd        = val;
   times->dom_bound_bck        = val;
   times->dom_posterior        = val;
   times->sp_biascorr          = val;
   times->dom_optacc           = val;
}

/* add times to running totals */
void 
WORK_times_add( WORKER* worker )
{
   TIMES* times = worker->times;
   TIMES* time_totals = worker->times_totals;

   /* totals */
   time_totals->program              += times->program;
   /* main loop */
   time_totals->loop                 += times->loop;
   /* indexes */
   time_totals->load_target_index    += times->load_target_index;
   time_totals->load_query_index     += times->load_query_index;
   /* load */
   time_totals->load_target          += times->load_target;
   time_totals->load_query           += times->load_query;
   /* naive algs */
   time_totals->naive_cloud          += times->naive_cloud;        
   /* quadratic algs */
   time_totals->quad_fwd             += times->quad_fwd;
   time_totals->quad_bck             += times->quad_bck;
   time_totals->quad_vit             += times->quad_vit;
   time_totals->quad_trace           += times->quad_trace;
   time_totals->quad_cloud_fwd       += times->quad_cloud_fwd;
   time_totals->quad_cloud_bck       += times->quad_cloud_bck;
   time_totals->quad_merge           += times->quad_merge;
   time_totals->quad_reorient        += times->quad_reorient;
   time_totals->quad_bound_fwd       += times->quad_bound_fwd;
   time_totals->quad_bound_bck       += times->quad_bound_bck;
   time_totals->quad_posterior       += times->quad_posterior;
   time_totals->quad_optacc          += times->quad_optacc;
   /* linear algs */   
   time_totals->lin_fwd              += times->lin_fwd;
   time_totals->lin_bck              += times->lin_bck;
   time_totals->lin_vit              += times->lin_vit;
   time_totals->lin_trace            += times->lin_trace;
   time_totals->lin_cloud_fwd        += times->lin_cloud_fwd;
   time_totals->lin_cloud_bck        += times->lin_cloud_bck ;
   time_totals->lin_merge            += times->lin_merge;
   time_totals->lin_reorient         += times->lin_reorient;
   time_totals->lin_bound_fwd        += times->lin_bound_fwd;
   time_totals->lin_bound_bck        += times->lin_bound_bck;
   /* sparse algs */
   time_totals->sp_build_mx          += times->sp_build_mx;
   time_totals->sp_fwd               += times->sp_fwd;
   time_totals->sp_bck               += times->sp_bck;
   time_totals->sp_vit               += times->sp_vit;
   time_totals->sp_trace             += times->sp_trace;
   time_totals->sp_cloud_fwd         += times->sp_cloud_fwd;
   time_totals->sp_cloud_bck         += times->sp_cloud_bck;
   time_totals->sp_bound_fwd         += times->sp_bound_fwd;
   time_totals->sp_bound_bck         += times->sp_bound_bck;
   time_totals->sp_posterior         += times->sp_posterior;
   time_totals->sp_biascorr          += times->sp_biascorr;
   time_totals->sp_optacc            += times->sp_optacc;
   /* doms */
   time_totals->dom_total            += times->dom_total;
   time_totals->dom_bound_fwd        += times->dom_bound_fwd;
   time_totals->dom_bound_bck        += times->dom_bound_bck;
   time_totals->dom_posterior        += times->dom_posterior;
   time_totals->dom_biascorr         += times->dom_biascorr;
   time_totals->dom_optacc           += times->dom_optacc;
}

