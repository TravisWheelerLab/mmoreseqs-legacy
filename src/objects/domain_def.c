/*******************************************************************************
 *  FILE:      dom_def.c
 *  PURPOSE:   DOMAIN_DEF object.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/* local imports */
#include "structs.h"
#include "../utilities/_utilities.h"
#include "_objects.h"

/* header */
#include "domain_def.h"

/*! FUNCTION:  DOMAIN_DEF_Create()
 *  SYNOPSIS:  Create new DOMAIN_DEF object and returns pointer.
 *             Most data is left NULL to be supplied by WORK_init().
 */
DOMAIN_DEF* 
DOMAIN_DEF_Create()
{
   DOMAIN_DEF* dom_def = NULL;
   dom_def = (DOMAIN_DEF*) ERROR_malloc( sizeof(DOMAIN_DEF) );

   /* domains vector */
   dom_def->N           = 0;
   dom_def->Nalloc      = 8;
   dom_def->domains     = ERROR_malloc( sizeof(DOMAIN_X) * dom_def->Nalloc );
   dom_def->dom_ranges  = VECTOR_RANGE_Create();
   dom_def->dom_fwdsc   = VECTOR_FLT_Create();
   dom_def->dom_bias    = VECTOR_FLT_Create();

   /* domain finding */
   dom_def->b_tot       = VECTOR_FLT_Create();
   dom_def->e_tot       = VECTOR_FLT_Create();
   dom_def->m_occ       = VECTOR_FLT_Create();

   /* null2 compo bias */
   dom_def->null2_sc    = VECTOR_FLT_Create();
   dom_def->null2_exp   = VECTOR_FLT_Create();
   dom_def->st_freq     = MATRIX_2D_Create(1, 1);
   dom_def->sp_freq     = VECTOR_FLT_Create();
   dom_def->st_num      = VECTOR_FLT_Create();

   /* alignment */
   dom_def->align       = ALIGNMENT_Create();
   dom_def->trace_aln   = VECTOR_TRACE_Create();
   dom_def->cigar_aln   = VECTOR_CHAR_Create();

   /* edgebounds */
   dom_def->edg         = EDGEBOUNDS_Create();

   /* default values */
   dom_def->rt1         = 0.25;
   dom_def->rt2         = 0.10;
   dom_def->rt3         = 0.20;

   /* totals */
   dom_def->dom_sumsc   = 0.0f;

   /* starting values */
   dom_def->n_regions      = 0;
   dom_def->n_domains      = 0;
   dom_def->n_envelopes    = 0;
   dom_def->n_clustered    = 0;
   dom_def->n_overlaps     = 0;

   return dom_def;
}

/*! FUNCTION:  DOMAIN_DEF_Reuse()
 *  SYNOPSIS:  Frees DOMAIN_DEF object and returns NULL pointer.
 */
DOMAIN_DEF* 
DOMAIN_DEF_Destroy( DOMAIN_DEF* dom_def )
{
   if (dom_def == NULL) return dom_def;

   /* domains vector */
   dom_def->domains     = ERROR_free( dom_def->domains );
   dom_def->dom_ranges  = VECTOR_RANGE_Destroy( dom_def->dom_ranges );
   dom_def->dom_fwdsc   = VECTOR_FLT_Destroy( dom_def->dom_fwdsc );
   dom_def->dom_bias    = VECTOR_FLT_Destroy( dom_def->dom_bias );

   /* domain finding */
   dom_def->b_tot       = VECTOR_FLT_Destroy( dom_def->b_tot );
   dom_def->e_tot       = VECTOR_FLT_Destroy( dom_def->e_tot );
   dom_def->m_occ       = VECTOR_FLT_Destroy( dom_def->m_occ );

   /* null2 compo bias */
   dom_def->null2_sc    = VECTOR_FLT_Destroy( dom_def->null2_sc );
   dom_def->null2_exp   = VECTOR_FLT_Destroy( dom_def->null2_exp );
   dom_def->st_freq     = MATRIX_2D_Destroy( dom_def->st_freq );
   dom_def->sp_freq     = VECTOR_FLT_Destroy( dom_def->sp_freq );
   dom_def->st_num      = VECTOR_FLT_Destroy( dom_def->st_num );

   /* alignment */
   dom_def->align       = ALIGNMENT_Destroy( dom_def->align );
   dom_def->trace_aln   = VECTOR_TRACE_Destroy( dom_def->trace_aln );
   dom_def->cigar_aln   = VECTOR_CHAR_Destroy( dom_def->cigar_aln );

   /* edgebounds */
   dom_def->edg         = EDGEBOUNDS_Destroy( dom_def->edg );

   dom_def = ERROR_free( dom_def );

   return dom_def;
}

/*! FUNCTION:  DOMAIN_DEF_Reuse()
 *  SYNOPSIS:  
 */
int 
DOMAIN_DEF_Reuse(   DOMAIN_DEF*    dom_def )
{
   /* domains vector */
   VECTOR_RANGE_Reuse( dom_def->dom_ranges );
   VECTOR_FLT_Reuse( dom_def->dom_fwdsc );
   VECTOR_FLT_Reuse( dom_def->dom_bias );

   /* domain finding */
   VECTOR_FLT_Reuse( dom_def->b_tot );
   VECTOR_FLT_Reuse( dom_def->e_tot );
   VECTOR_FLT_Reuse( dom_def->m_occ );

   /* null2 compo bias */
   VECTOR_FLT_Reuse( dom_def->null2_sc );
   VECTOR_FLT_Reuse( dom_def->null2_exp );
   MATRIX_2D_Reuse( dom_def->st_freq, 1, 1 );
   VECTOR_FLT_Reuse( dom_def->sp_freq );
   VECTOR_FLT_Reuse( dom_def->st_num );

   /* alignment */
   ALIGNMENT_Reuse( dom_def->align, 0, 0 );
   VECTOR_TRACE_Reuse( dom_def->trace_aln );
   VECTOR_CHAR_Reuse( dom_def->cigar_aln );

   return STATUS_SUCCESS;
}

/*! FUNCTION:  DOMAIN_DEF_Resize()
 *  SYNOPSIS:  
 */
int
DOMAIN_DEF_Resize(   DOMAIN_DEF*    dom_def,
                     int            size )
{
   VECTOR_FLT_GrowTo( dom_def->b_tot, size );
   VECTOR_FLT_GrowTo( dom_def->e_tot, size );
   VECTOR_FLT_GrowTo( dom_def->m_occ, size );
   VECTOR_FLT_GrowTo( dom_def->null2_sc, size );

   return STATUS_SUCCESS;
}

/*! FUNCTION:  DOMAIN_DEF_GrowTo()
 *  SYNOPSIS:  
 */
int
DOMAIN_DEF_GrowTo(   DOMAIN_DEF*    dom_def,
                     int            size )
{
   if ( dom_def->b_tot->N < size ) {
      DOMAIN_DEF_Resize( dom_def, size );
   }
   return STATUS_SUCCESS;
}


/*! FUNCTION:  DOMAIN_DEF_Pushback()
 *  SYNOPSIS:  
 */
int
DOMAIN_DEF_Pushback(   DOMAIN_DEF*    dom_def,
                       DOMAIN_X       domain )
{
   dom_def->N += 1;
   if ( dom_def->N >= dom_def->Nalloc ) {
      DOMAIN_DEF_Resize( dom_def, dom_def->N * 2 );
   }
   dom_def->domains[dom_def->N] = domain;
}