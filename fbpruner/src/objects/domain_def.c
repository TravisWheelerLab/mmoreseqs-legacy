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
#include "../utilities/utilities.h"
#include "objects.h"

/* header */
#include "domain_def.h"

/** FUNCTION:  DOMAIN_DEF_Create()
 *  SYNOPSIS:  Create new DOMAIN_DEF object and returns pointer.
 *             Most data is left NULL to be supplied by WORK_init().
 */
DOMAIN_DEF* 
DOMAIN_DEF_Create()
{
   DOMAIN_DEF* dom_def = NULL;
   
   dom_def = (DOMAIN_DEF*) ERROR_malloc( sizeof(DOMAIN_DEF) );

   dom_def->b_tot    = VECTOR_FLT_Create();
   dom_def->e_tot    = VECTOR_FLT_Create();
   dom_def->m_occ    = VECTOR_FLT_Create();
   dom_def->null2_sc = VECTOR_FLT_Create();
   dom_def->st_freq  = MATRIX_2D_Create(1, 1);
   dom_def->sp_freq  = VECTOR_FLT_Create();

   /* default values */
   dom_def->rt1           = 0.25;
   dom_def->rt2           = 0.10;
   dom_def->rt3           = 0.20;

   /* empty values */

   return dom_def;
}

/** FUNCTION:  DOMAIN_DEF_Destroy()
 *  SYNOPSIS:  Frees DOMAIN_DEF object and returns pointer.
 */
DOMAIN_DEF* 
DOMAIN_DEF_Destroy( DOMAIN_DEF* dom_def )
{
   if (dom_def == NULL) return dom_def;

   VECTOR_FLT_Destroy( dom_def->b_tot );
   VECTOR_FLT_Destroy( dom_def->e_tot );
   VECTOR_FLT_Destroy( dom_def->m_occ );
   VECTOR_FLT_Destroy( dom_def->null2_sc );
   MATRIX_2D_Destroy( dom_def->st_freq );
   VECTOR_FLT_Destroy( dom_def->sp_freq );

   dom_def = ERROR_free( dom_def );

   return dom_def;
}

/** FUNCTION:  DOMAIN_DEF_Reuse()
 *  SYNOPSIS:  
 */
void
DOMAIN_DEF_Reuse(   DOMAIN_DEF*    dom_def )
{
   VECTOR_FLT_Reuse( dom_def->b_tot );
   VECTOR_FLT_Reuse( dom_def->e_tot );
   VECTOR_FLT_Reuse( dom_def->m_occ );
   VECTOR_FLT_Reuse( dom_def->null2_sc );
}

/** FUNCTION:  DOMAIN_DEF_Resize()
 *  SYNOPSIS:  
 */
void
DOMAIN_DEF_Resize(   DOMAIN_DEF*    dom_def,
                     int            size )
{
   VECTOR_FLT_GrowTo( dom_def->b_tot, size );
   VECTOR_FLT_GrowTo( dom_def->e_tot, size );
   VECTOR_FLT_GrowTo( dom_def->m_occ, size );
   VECTOR_FLT_GrowTo( dom_def->null2_sc, size );
}

/** FUNCTION:  DOMAIN_DEF_GrowTo()
 *  SYNOPSIS:  
 */
void
DOMAIN_DEF_GrowTo(   DOMAIN_DEF*    dom_def,
                     int            size )
{
   if ( dom_def->b_tot->N < size ) {
      DOMAIN_DEF_Resize( dom_def, size );
   }
}

