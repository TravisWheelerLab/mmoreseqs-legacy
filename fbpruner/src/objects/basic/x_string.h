/*******************************************************************************
 *  FILE:      x_string.h
 *  PURPOSE:   X_STRING object
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _X_STRING_H
#define _X_STRING_H

/*
 *  FUNCTION:  X_STRING_Create()
 *  SYNOPSIS:
 */
 X_STRING* 
 X_STRING_Create( char* chars );

/*
 *  FUNCTION:  X_STRING_Destroy()
 *  SYNOPSIS:
 */
X_STRING* 
X_STRING_Destroy( X_STRING* str );

#endif /* _X_STRING_H */