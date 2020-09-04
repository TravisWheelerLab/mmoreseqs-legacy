/*******************************************************************************
 *  FILE:      debugger.c
 *  PURPOSE:   DEBUGGER Object.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:     
 *******************************************************************************/

#ifndef _DEBUGGER_H
#define _DEBUGGER_H

/*
 *  FUNCTION:  DEBUGGER_Create()
 *  SYNOPSIS:
 */
DEBUG_KIT* DEBUGGER_Create();

/*
 *  FUNCTION:  DEBUGGER_Destroy()
 *  SYNOPSIS:
 */
void* DEBUGGER_Destroy( DEBUG_KIT* dbg );

#endif /* _DEBUGGER_H */