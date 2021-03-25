/*******************************************************************************
 *  FILE:      debugger.c
 *  PURPOSE:   DEBUGGER Object.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:     
 *******************************************************************************/

#ifndef _DEBUGGER_H
#define _DEBUGGER_H

/*  FUNCTION:     DEBUGGER_Create()
 *  SYNOPSIS:
 */
DEBUG_KIT* DEBUGGER_Create( char*   filepath );

/*  FUNCTION:     DEBUGGER_Destroy()
 *  SYNOPSIS:     
 */
void* DEBUGGER_Destroy( DEBUG_KIT*  dbg );

/*  FUNCTION:     DEBUGGER_Reuse()
 *  SYNOPSIS:     
 */
void* DEBUGGER_Reuse(   DEBUG_KIT*  dbg,
                        int         Q,
                        int         T );

/*  FUNCTION:     DEBUGGER_Make_Dir()
 *  SYNOPSIS:     Attempt to create debugger folder.
 */
 void DEBUGGER_Make_Dir(   DEBUG_KIT*  dbg );

#endif /* _DEBUGGER_H */