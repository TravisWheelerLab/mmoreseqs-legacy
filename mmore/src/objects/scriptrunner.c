/*******************************************************************************
 *  FILE:      scriptrunner.c
 *  PURPOSE:   SCRIPTRUNNER object.
 *             Sets up and executes scripts.
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

/* header */
#include "_objects.h"
#include "scriptrunner.h"

/*! FUNCTION:  SCRIPTRUNNER_Create()
 *  SYNOPSIS:  Create new SCRIPTRUNNER object, allocates memory and returns pointer.
 */
SCRIPTRUNNER* 
SCRIPTRUNNER_Create()
{
   SCRIPTRUNNER* runner;
   runner = ERROR_malloc( sizeof(SCRIPTRUNNER) );

   runner->tool         = NULL;
   runner->script_path  = NULL;
   runner->env_names    = VECTOR_STR_Create();
   runner->env_values   = VECTOR_STR_Create();
   runner->arg_flags    = VECTOR_STR_Create();
   runner->arg_values   = VECTOR_STR_Create();
   runner->command      = VECTOR_STR_Create();

   return runner;
}

/*! FUNCTION:  SCRIPTRUNNER_Destroy()
 *  SYNOPSIS:  Deletes SCRIPTRUNNER and frees memory.
 */
SCRIPTRUNNER* 
SCRIPTRUNNER_Destroy( SCRIPTRUNNER*    runner )
{
   if ( runner == NULL ) return runner;

   runner->tool         = STR_Destroy( runner->tool );
   runner->script_path  = STR_Destroy( runner->script_path );
   runner->env_names    = VECTOR_STR_Destroy( runner->env_names );
   runner->env_values   = VECTOR_STR_Destroy( runner->env_values );
   runner->arg_flags    = VECTOR_STR_Destroy( runner->arg_flags );
   runner->arg_values   = VECTOR_STR_Destroy( runner->arg_values );
   runner->command      = VECTOR_STR_Destroy( runner->command );

   runner               = ERROR_free( runner );

   return runner;
}

/*! FUNCTION:  SCRIPTRUNNER_SetTool()
 *  SYNOPSIS:  Sets which tool (bash, python, etc.) that <runner> will use.
 *             If <script> is an executable, can be left NULL;
 */
STATUS_FLAG 
SCRIPTRUNNER_SetTool(  SCRIPTRUNNER*     runner,
                        STR               tool )
{
   runner->tool = STR_Create( tool );

   return STATUS_SUCCESS;
}

/*! FUNCTION:  SCRIPTRUNNER_SetFilepath()
 *  SYNOPSIS:  Sets the filepath of the script to be run by <runner>.
 */
STATUS_FLAG 
SCRIPTRUNNER_SetScript(    SCRIPTRUNNER*     runner,
                           STR               script_path )
{
   runner->script_path = STR_Create(script_path);

   return STATUS_SUCCESS;
}

/*! FUNCTION:  SCRIPTRUNNER_Add_Env_Variable()
 *  SYNOPSIS:  Adds environmental variable to <runner>.
 */
STATUS_FLAG 
SCRIPTRUNNER_Add_Env_Variable(   SCRIPTRUNNER*     runner,
                                 STR               env_name,
                                 STR               env_value )
{
   /* if value is NULL, does not create a variable. */
   if ( env_name == NULL ) {
      fprintf(stderr, "ERROR: env_name is NULL. ( env_name = %s, env_value = %s )\n", env_name, env_value );
   }
   if ( env_value == NULL ) {
      return STATUS_FAILURE;
   } 
   
   VECTOR_STR_Pushback( runner->env_names, env_name );
   VECTOR_STR_Pushback( runner->env_values, env_value );

   return STATUS_SUCCESS;
}

/*! FUNCTION:  SCRIPTRUNNER_If_Add_Env_Variable()
 *  SYNOPSIS:  Adds environmental variable to <runner> if condition is true.
 */
STATUS_FLAG 
SCRIPTRUNNER_If_Add_Env_Variable(   SCRIPTRUNNER*     runner,
                                    STR               env_name,
                                    STR               env_value,
                                    bool              condition )
{
   /* if name or value is NULL, does not create a variable. */
   if ( env_name == NULL || env_value == NULL ) {
      return STATUS_FAILURE;
   } 
   
   if ( condition == true ) {
      VECTOR_STR_Pushback( runner->env_names, env_name );
      VECTOR_STR_Pushback( runner->env_values, env_value );
   }

   return STATUS_SUCCESS;
}

/*! FUNCTION:  SCRIPTRUNNER_Add_Script_Argument()
 *  SYNOPSIS:  Adds script argument to <runner>.
 */
STATUS_FLAG 
SCRIPTRUNNER_Add_Script_Argument(   SCRIPTRUNNER*     runner,
                                    STR               arg_flag,
                                    STR               arg_value )
{
   /* if value is null, do not add argument */
   if ( arg_value != NULL ) {
      VECTOR_STR_Pushback( runner->arg_flags, arg_flag );
      VECTOR_STR_Pushback( runner->arg_values, arg_value );

      return STATUS_SUCCESS;
   }

   return STATUS_FAILURE;
}

/*! FUNCTION:  SCRIPTRUNNER_Check()
 *  SYNOPSIS:  Verifies that <runner>  args are valid.
 */
STATUS_FLAG 
SCRIPTRUNNER_Check(    SCRIPTRUNNER*     runner )
{
   /* check if tool is installed on system */

   /* check if script path exists */
   bool file_exists = SYSTEMIO_FileExists( runner->script_path );
   if ( file_exists == false ) {
      fprintf( stderr, "ERROR: script '%s' for SCRIPTRUNNER does not exist.\n", 
         runner->script_path );
      return STATUS_FAILURE;
   }
   return STATUS_SUCCESS;
}

/*! FUNCTION:  SCRIPTRUNNER_Execute()
 *  SYNOPSIS:  Run <runner> script with given arguments.
 */
STATUS_FLAG 
SCRIPTRUNNER_Execute(    SCRIPTRUNNER*     runner )
{
   int num_vars;

   VECTOR_STR_Reuse( runner->command );

   /* check if arguments are valid */
   SCRIPTRUNNER_Check( runner );

   /* load environmental variables */
   num_vars = VECTOR_STR_GetSize( runner->env_names ); 
   for (int i = 0; i < num_vars; i++ )  
   {
      STR env_name   = VEC_X( runner->env_names, i );
      STR env_val    = VEC_X( runner->env_values, i );
      SYSTEMIO_AddEnvironmentalVar( env_name, env_val );
   }

   /* load tool and filepath */
   if ( runner->tool != NULL ) {
      VECTOR_STR_Pushback( runner->command, runner->tool );
   }
   /* load filepath to script */
   VECTOR_STR_Pushback( runner->command, runner->script_path );

   /* load arguments into command */
   num_vars = VECTOR_STR_GetSize( runner->arg_flags );
   for (int i = 0; i < num_vars; i++ ) 
   {
      STR arg_flag   = VEC_X( runner->arg_flags, i );
      STR arg_val    = VEC_X( runner->arg_values, i );
      if ( arg_flag != NULL ) {
         VECTOR_STR_Pushback( runner->command, arg_flag );
      }
      VECTOR_STR_Pushback( runner->command, arg_val );
   }
   /* array must end with NULL pointer */
   VECTOR_STR_Pushback( runner->command, NULL );

   /* extract array from vector */
   // VECTOR_STR_Dump( runner->command, stdout );
   STR* command_array = VECTOR_STR_GetArray( runner->command );
   int ERRORCHECK_exit_code = execvp( command_array[0], command_array );

   /* program should not get here */
   fprintf( stderr, "ERROR: SCRIPTRUNNER failed with code (%d).\n", ERRORCHECK_exit_code );

   return STATUS_SUCCESS;
}