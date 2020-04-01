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

/* destructor */
void WORKER_Destory( WORKER* work );

#endif /* _WORKER_H */