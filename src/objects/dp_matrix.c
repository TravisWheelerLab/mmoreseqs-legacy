/*******************************************************************************
 *  - FILE:      dp_mx.c
 *  - DESC:    DP_MATRIX object.
 *             Master dp_mx object. Maintains memory for all data structures used in pipeline.
 *             Contains data shared by all dp_mx threads.
 *  NOTES:
 *    - WIP.
 *    - This DP_MATRIX is designed to support the different types of matrices, determined by DPMX_MODE.
 *    - This should use a mutually-shared databank, so the data can be used for
 *      any DP_MATRIX needed.
 *    - Will require a Create() function for each that does not allocate matrix data.
 *    - Then Nalloc and ->data location must be passed between matrices.
 *    - Same with Destroy() function, watch out for double-free.
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
#include "dp_matrix.h"

/** FUNCTION:  DP_MATRIX_Create()
 *  SYNOPSIS:  Create new DP_MATRIX object and returns pointer.
 *             Most data is left NULL to be supplied by WORK_init().
 */
DP_MATRIX*
DP_MATRIX_Create(DPMX_MODE mode) {
  DP_MATRIX* dp_mx;
  dp_mx = ERROR_malloc(sizeof(DP_MATRIX));

  /* matrix types */
  dp_mx->mode = mode;

  /* Set unneccessary matrices to null */
  dp_mx->st_MX = NULL;
  dp_mx->st_MX3 = NULL;
  dp_mx->st_SMX = NULL;
  dp_mx->sp_MX = NULL;

  /* build necessary matrices */
  dp_mx->st_MX = MATRIX_3D_Create(1, 1, 1);
  dp_mx->st_MX3 = MATRIX_3D_Create(1, 1, 1);
  dp_mx->st_SMX = MATRIX_3D_SPARSE_Create();

  dp_mx->sp_MX = MATRIX_2D_Create(1, 1);

  /* default values */
  dp_mx->Q = 0;
  dp_mx->T = 0;

  return dp_mx;
}

/** FUNCTION:  DP_MATRIX_Reuse()
 *  SYNOPSIS:
 */
DP_MATRIX*
DP_MATRIX_Reuse(DP_MATRIX* dp_mx,
                int Q,
                int T,
                DPMX_MODE mode) {
}

/** FUNCTION:  DP_MATRIX_Destroy()
 *  SYNOPSIS:  Frees DP_MATRIX object and returns pointer.
 */
DP_MATRIX*
DP_MATRIX_Destroy(DP_MATRIX* dp_mx) {
  if (dp_mx == NULL)
    return dp_mx;

  dp_mx->st_MX = MATRIX_3D_Destroy(dp_mx->st_MX);
  dp_mx->st_MX3 = MATRIX_3D_Destroy(dp_mx->st_MX3);
  dp_mx->st_SMX = MATRIX_3D_SPARSE_Destroy(dp_mx->st_SMX);
  dp_mx->sp_MX = MATRIX_2D_Destroy(dp_mx->sp_MX);

  dp_mx = ERROR_free(dp_mx);

  return dp_mx;
}
