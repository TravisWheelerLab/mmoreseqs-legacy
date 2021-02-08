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
   runner               = ERROR_free( runner );

   return runner;
}

/*! FUNCTION:  SCRIPTRUNNER_SetTool()
 *  SYNOPSIS:  Sets which tool (bash, python, etc.) that <runner> will use.
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
SCRIPTRUNNER_SetScript(   SCRIPTRUNNER*     runner,
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
   VECTOR_STR_Pushback( runner->env_names, env_name );
   VECTOR_STR_Pushback( runner->env_values, env_value );

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
   VECTOR_STR_Pushback( runner->arg_flags, arg_flag );
   VECTOR_STR_Pushback( runner->arg_values, arg_value );

   return STATUS_SUCCESS;
}

/*! FUNCTION:  SCRIPTRUNNER_Execute()
 *  SYNOPSIS:  Run <runner> script with given arguments.
 */
STATUS_FLAG 
SCRIPTRUNNER_Execute(    SCRIPTRUNNER*     runner )
{
   
}