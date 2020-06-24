/*******************************************************************************
 *  FILE:      worker.c
 *  PURPOSE:   WORKER object.h
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _WORKER_H
#define _WORKER_H

/* constructor */
WORKER* WORKER_Create();

/* constructor with args supplied */
WORKER* WORKER_Create_with_Args( ARGS* args );

/* destructor */
void* WORKER_Destroy( WORKER* worker );

#endif /* _WORKER_H */