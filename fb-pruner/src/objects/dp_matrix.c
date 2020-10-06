/*******************************************************************************
 *  FILE:      dp_mx.c
 *  PURPOSE:   DP_MATRIX object.
 *             Master dp_mx object. Maintains memory for all data structures used in pipeline.
 *             Contains data shared by all dp_mx threads.
 *             TODO: Worker will have thread dp_mxs.
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
#include "utilities.h"
#include "objects.h"

/* header */
#include "dp_matrix.h"

/** FUNCTION:  DP_MATRIX_Create()
 *  SYNOPSIS:  Create new DP_MATRIX object and returns pointer.
 *             Most data is left NULL to be supplied by WORK_init().
 */
DP_MATRIX* 
DP_MATRIX_Create( int   mx_type )
{
   DP_MATRIX* dp_mx;
   
   dp_mx = ERROR_malloc( sizeof(DP_MATRIX) );

   /* build necessary matrices */
   dp_mx->st_MX   = NULL;
   dp_mx->st_MX3  = NULL;
   dp_mx->st_SMX  = NULL;
   dp_mx->sp_MX   = NULL;

   /* build necessary matrices */
   dp_mx->st_MX   = MATRIX_3D_Create(1, 1, 1);
   dp_mx->st_MX3  = MATRIX_3D_Create(1, 1, 1);
   dp_mx->st_SMX  = MATRIX_3D_SPARSE_Create();

   dp_mx->sp_MX   = MATRIX_2D_Create(1, 1);

   /* default values */
   dp_mx->Q       = 0;
   dp_mx->T       = 0;
   dp_mx->mx_type = false;

   return dp_mx;
}

/** FUNCTION:  DP_MATRIX_GrowTo()
 *  SYNOPSIS:  
 */
DP_MATRIX* 
DP_MATRIX_GrowTo(    DP_MATRIX*     dp_mx,
                     int            size )
{

}

/** FUNCTION:  DP_MATRIX_Destroy()
 *  SYNOPSIS:  Frees DP_MATRIX object and returns pointer.
 */
DP_MATRIX* 
DP_MATRIX_Destroy( DP_MATRIX* dp_mx )
{
   if (dp_mx == NULL) return dp_mx;
   
   dp_mx->st_MX   = MATRIX_3D_Destroy( dp_mx->st_MX );
   dp_mx->st_MX3  = MATRIX_3D_Destroy( dp_mx->st_MX3 );
   dp_mx->st_SMX  = MATRIX_3D_SPARSE_Destroy( dp_mx->st_SMX );
   dp_mx->sp_MX   = MATRIX_2D_Destroy( dp_mx->sp_MX );

   dp_mx = ERROR_free( dp_mx );

   return dp_mx;
}