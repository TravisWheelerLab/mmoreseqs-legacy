/*******************************************************************************
 *  - FILE:      scriptrunner.c
 *  - DESC:    SCRIPTRUNNER object.
 *             Sets up and executes scripts.
 *******************************************************************************/

#ifndef _SCRIPTRUNNER_H
#define _SCRIPTRUNNER_H

/*! FUNCTION:  SCRIPTRUNNER_Create()
 *  SYNOPSIS:  Create new SCRIPTRUNNER object, allocates memory and returns
 * pointer.
 */
SCRIPTRUNNER* SCRIPTRUNNER_Create();

/*! FUNCTION:  SCRIPTRUNNER_Destroy()
 *  SYNOPSIS:  Deletes SCRIPTRUNNER and frees memory.
 */
SCRIPTRUNNER* SCRIPTRUNNER_Destroy(SCRIPTRUNNER* runner);

/*! FUNCTION:  SCRIPTRUNNER_SetTool()
 *  SYNOPSIS:  Sets which tool (bash, python, etc.) that <runner> will use.
 */
STATUS_FLAG
SCRIPTRUNNER_SetTool(SCRIPTRUNNER* runner, STR tool);
/*! FUNCTION:  SCRIPTRUNNER_SetFilepath()
 *  SYNOPSIS:  Sets the filepath of the script to be run by <runner>.
 */
STATUS_FLAG
SCRIPTRUNNER_SetScript(SCRIPTRUNNER* runner, STR script_path);

/*! FUNCTION:  SCRIPTRUNNER_Add_Env_Variable()
 *  SYNOPSIS:  Adds environmental variable to <runner>.
 */
STATUS_FLAG
SCRIPTRUNNER_Add_Env_Variable(SCRIPTRUNNER* runner,
                              STR env_name,
                              STR env_value);

/*! FUNCTION:  SCRIPTRUNNER_If_Add_Env_Variable()
 *  SYNOPSIS:  Adds environmental variable to <runner> if condition is true.
 */
STATUS_FLAG
SCRIPTRUNNER_If_Add_Env_Variable(SCRIPTRUNNER* runner,
                                 STR env_name,
                                 STR env_value,
                                 bool condition);

/*! FUNCTION:  SCRIPTRUNNER_Add_Script_Argument()
 *  SYNOPSIS:  Adds script argument to <runner>.
 */
STATUS_FLAG
SCRIPTRUNNER_Add_Script_Argument(SCRIPTRUNNER* runner,
                                 STR arg_flag,
                                 STR arg_value);

/*! FUNCTION:  SCRIPTRUNNER_Execute()
 *  SYNOPSIS:  Run <runner> script with given arguments.
 */
STATUS_FLAG
SCRIPTRUNNER_Execute(SCRIPTRUNNER* runner);

/*! FUNCTION:  SCRIPTRUNNER_Dump()
 *  SYNOPSIS:  Print out script data.
 */
STATUS_FLAG
SCRIPTRUNNER_Dump(SCRIPTRUNNER* runner, FILE* fp);

#endif /* _SCRIPTRUNNER_H */
