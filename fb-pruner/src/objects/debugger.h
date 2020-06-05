/*******************************************************************************
 *  FILE:      debugger.c
 *  PURPOSE:   DEBUGGER Object.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:     
 *******************************************************************************/

#ifndef _DEBUGGER_H
#define _DEBUGGER_H

/* constructor */
DEBUG_KIT* DEBUGGER_Create();

/* destructor */
void DEBUGGER_Destroy( DEBUG_KIT* debugger );

#endif /* _DEBUGGER_H */