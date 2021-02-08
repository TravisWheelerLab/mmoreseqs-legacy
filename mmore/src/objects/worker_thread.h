/*******************************************************************************
 *  FILE:   	worker_thread.h
 *  PURPOSE:   WORKER_THREAD object
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _WORKER_THREAD_H
#define _WORKER_THREAD_H

/* constructor */
WORKER_THREAD* 
WORKER_THREAD_Create();

/* destructor */
void* 
WORKER_THREAD_Destroy( WORKER_THREAD* worker );

#endif /* _WORKER_H */