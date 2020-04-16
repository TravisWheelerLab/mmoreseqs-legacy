/*******************************************************************************
 *  FILE:      edgebound.c
 *  PURPOSE:   EDGEBOUNDS Object
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _EDGEBOUND_H
#define _EDGEBOUND_H

/* constructor */
EDGEBOUNDS* EDGEBOUNDS_Create();

/* initializer */
void EDGEBOUNDS_Init(EDGEBOUNDS** edg);

/* destructor */
void EDGEBOUNDS_Destroy(EDGEBOUNDS* edg);

/* reuse (sets length to 0 by ) */
void EDGEBOUNDS_Reuse(EDGEBOUNDS*  edg);

/* access bound in edgebound array */
BOUND* EDGEBOUNDS_Get(EDGEBOUNDS*   edg, 
                      int           i );

/* push BOUND onto EDGEBOUNDS array */
void EDGEBOUNDS_Pushback(EDGEBOUNDS* edg,
                         BOUND       bnd);

/* push HEAD onto EDGEBOUNDS array */
void EDGEBOUNDS_Pushback_Head(EDGEBOUNDS* edg,
                              int         id,
                              int         head);

/* insert BOUND into index of EDGEBOUNDS array */
void EDGEBOUNDS_Insert(EDGEBOUNDS* edg,
                       BOUND       bnd,
                       int         i);

/* insert BOUND into delete of EDGEBOUNDS array */
void EDGEBOUNDS_Delete(EDGEBOUNDS* edg,
                       BOUND       bnd,
                       int         i);

/* remove all BOUNDS from EDGEBOUND */
void EDGEBOUNDS_Clear(EDGEBOUNDS* edg);

/* resize BOUNDS array */
void EDGEBOUNDS_Resize(EDGEBOUNDS* edg, 
                       int         size);

/* reverse ordering of BOUNDS array */
void EDGEBOUNDS_Reverse(EDGEBOUNDS* edg);

/* find the HEADS of each set of BOUNDS arrays by ID */
void EDGEBOUNDS_SetHeads(EDGEBOUNDS* edg);

/* Output EDGEBOUNDS to FILE POINTER */
void EDGEBOUNDS_Dump(EDGEBOUNDS*  edg,
                     FILE*        fp);

/* Output EDGEBOUNDS to FILE based on FILENAME */
void EDGEBOUNDS_Save(EDGEBOUNDS* edg,
                     const char* _filename_);

/* Sort BOUNDS array in EDGEBOUNDS */
void EDGEBOUNDS_Sort(EDGEBOUNDS *edg);

/* Count all cells covered by BOUNDS array in EDGEBOUNDS */
int EDGEBOUNDS_Count(EDGEBOUNDS *edg);

#endif /* _EDGEBOUND_H */