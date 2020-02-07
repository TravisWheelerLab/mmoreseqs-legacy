/*******************************************************************************
 *  @file edgebound.h
 *  @brief Edgebounds Datatype
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

#ifndef _EDGEBOUND_H
#define _EDGEBOUND_H

EDGEBOUNDS *edgebounds_Create();

void edgebounds_Init(EDGEBOUNDS **edg);

void edgebounds_Destroy(EDGEBOUNDS *edg);

void edgebounds_Add(EDGEBOUNDS *edg,
                    BOUND bnd);

void edgebounds_Insert(EDGEBOUNDS *edg,
                       BOUND bnd,
                       int i);

void edgebounds_Delete(EDGEBOUNDS *edg,
                       BOUND bnd,
                       int i);

void edgebounds_Resize(EDGEBOUNDS *edg);

void edgebounds_Reverse(EDGEBOUNDS *edg);

void edgebounds_Print(EDGEBOUNDS *edg);
void bound_Print(BOUND bnd);

void edgebounds_Save(EDGEBOUNDS *edg,
                     const char *_filename_);

void edgbounds_Sort(EDGEBOUNDS *edg);

int edgebounds_Count(EDGEBOUNDS *edg);

#endif /* _EDGEBOUND_H */