/*******************************************************************************
 *  FILE:      debugger.c
 *  PURPOSE:   DEBUGGER Object.
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
#include <ctype.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

/* local imports */
#include "structs.h"
#include "../utilities/utilities.h"
#include "objects.h"

/* header */
#include "debugger.h"

/*  FUNCTION:     DEBUGGER_Create()
 *  SYNOPSIS:
 */
DEBUG_KIT* DEBUGGER_Create( char*   filepath )
{
	DEBUG_KIT* dbg = NULL;
   dbg = (DEBUG_KIT*) ERROR_malloc( sizeof(DEBUG_KIT) );

   dbg->is_debugging 	   = false;
   dbg->verbose_level 	   = 0;
   dbg->is_embed 			   = false;
   dbg->is_viz 			   = false;

   dbg->dbgout_dir 	      = strdup( filepath );
   dbg->dbgout_fp 			= NULL;

   dbg->cloud_MX           = MATRIX_2D_Create( 1, 1 ); 
   dbg->cloud_MX3          = MATRIX_2D_Create( 1, 1 ); 

   dbg->test_MX            = MATRIX_3D_Create( NUM_NORMAL_STATES, 1, 1 );
   dbg->test_MX3           = MATRIX_3D_Create( NUM_NORMAL_STATES, 1, 1 );

   dbg->test_edg           = EDGEBOUNDS_Create();

   return dbg;
}

/*  FUNCTION:     DEBUGGER_Destroy()
 *  SYNOPSIS:     
 */
void* DEBUGGER_Destroy( DEBUG_KIT*  dbg )
{
	dbg->dbgout_dir      = ERROR_free( dbg->dbgout_dir );
	dbg->dbgout_fp  	   = ERROR_free( dbg->dbgout_fp );

   dbg->cloud_MX        = MATRIX_2D_Destroy( dbg->cloud_MX );
   dbg->cloud_MX3       = MATRIX_2D_Destroy( dbg->cloud_MX3 );

	dbg->test_MX 	      = MATRIX_3D_Destroy( dbg->test_MX );
	dbg->test_MX3 	      = MATRIX_3D_Destroy( dbg->test_MX3 );

	dbg->test_edg 	      = EDGEBOUNDS_Destroy( dbg->test_edg );

   dbg = ERROR_free( dbg );

   return NULL;
}

/*  FUNCTION:     DEBUGGER_Reuse()
 *  SYNOPSIS:     
 */
void* DEBUGGER_Reuse(   DEBUG_KIT*  dbg,
                        int         Q,
                        int         T )
{
   MATRIX_2D_Reuse( dbg->cloud_MX, Q+1, T+1 );
   MATRIX_2D_Reuse( dbg->cloud_MX3, 3, (Q+1)+(T+1) );

   MATRIX_3D_Reuse( dbg->test_MX, NUM_NORMAL_STATES, Q+1, T+1 );
   MATRIX_3D_Reuse( dbg->test_MX3, NUM_NORMAL_STATES, 3, (Q+1)+(T+1) );

   EDGEBOUNDS_Reuse( dbg->test_edg, Q, T );
}

/*  FUNCTION:     DEBUGGER_Make_Dir()
 *  SYNOPSIS:     Attempt to create debugger folder.
 */
 void DEBUGGER_Make_Dir(   DEBUG_KIT*  dbg ) 
 {
   int status = mkdir(dbg->dbgout_dir, 0777);
   if ( status != -1 ) {
      fprintf( stderr, "# Folder created successfully: %s\n", dbg->dbgout_dir );
   } else {
      fprintf( stderr, "# Folder already exists: %s\n", dbg->dbgout_dir );
   }
 }

