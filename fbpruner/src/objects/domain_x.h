/*******************************************************************************
 *  FILE:      domain_x.c
 *  PURPOSE:   DOMAIN_X object.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _DOMAIN_X_H
#define _DOMAIN_X_H

/** FUNCTION:  DOMAIN_X_Create()
 *  SYNOPSIS:  Create new DOMAIN_X object and returns pointer.
 *             Most data is left NULL to be supplied by WORK_init().
 */
DOMAIN_X* 
DOMAIN_X_Create();


/** FUNCTION:  DOMAIN_X_Destroy()
 *  SYNOPSIS:  Frees DOMAIN_X object and returns pointer.
 */
DOMAIN_X* 
DOMAIN_X_Destroy( DOMAIN_X* dom );


/** FUNCTION:  DOMAIN_X_Reuse()
 *  SYNOPSIS:  
 */
int
DOMAIN_X_Reuse(   DOMAIN_X*    dom );


#endif /* _DOMAIN_X_H */