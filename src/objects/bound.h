/*******************************************************************************
 *  FILE:      edgebound.c
 *  PURPOSE:   EDGEBOUNDS Object
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/


#ifndef _BOUND_H
#define _BOUND_H

void BOUND_Dump(const BOUND  bnd,
                FILE*        fp);

int BOUND_Compare(const BOUND  a, 
                  const BOUND  b);

#endif /* _BOUND_H */