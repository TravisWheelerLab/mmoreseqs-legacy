/*******************************************************************************
 *  FILE:      worker.c
 *  PURPOSE:   WORKER object.h
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _WORKER_H
#define _WORKER_H

/*
 *  FUNCTION:  WORKER_Create()
 *  SYNOPSIS:  Create new WORKER object and returns pointer.
 *             Most data is left NULL to be supplied by WORK_init().
 */
WORKER* WORKER_Create();

/*
 *  FUNCTION:  WORKER_Create_with_Args()
 *  SYNOPSIS:  Create new WORKER object and adds <args> to worker.
 */
WORKER* WORKER_Create_with_Args( ARGS* args );

/*
 *  FUNCTION:  WORKER_Create_Threads()
 *  SYNOPSIS:  Creates <N_threads> WORKER_THREAD objects for <worker>.
 */
void
WORKER_Create_Threads(  WORKER*  worker,
                        int      N_threads );

/*
 *  FUNCTION:  WORKER_Destroy()
 *  SYNOPSIS:  Frees WORKER object and returns pointer.
 */
void* 
WORKER_Destroy( WORKER* worker );

#endif /* _WORKER_H */