/*******************************************************************************
 *  FILE:      edgebound.c
 *  PURPOSE:   EDGEBOUNDS Object
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _EDGEBOUND_H
#define _EDGEBOUND_H

/*
 *  FUNCTION:  EDGEBOUNDS_Create()
 *  SYNOPSIS:  Create new EDGEBOUNDS object and returns pointer.
 */
EDGEBOUNDS* EDGEBOUNDS_Create( void );

/*
 *  FUNCTION:  EDGEBOUNDS_Create()
 *  SYNOPSIS:  Create new EDGEBOUNDS object with chosen size and returns pointer.
 */
EDGEBOUNDS* EDGEBOUNDS_Create_by_Size( const int size );

/*
 *  FUNCTION: EDGEBOUNDS_Destroy()
 *  SYNOPSIS: Frees all memory from EDGEBOUNDS object.
 */
void* EDGEBOUNDS_Destroy( EDGEBOUNDS*  edg );

/*
 *  FUNCTION: EDGEBOUNDS_Reuse()
 *  SYNOPSIS: Reuses EDGEBOUNDS by "clearing" edgebound list (does not realloc).
 */
void EDGEBOUNDS_Reuse( EDGEBOUNDS*   edg, 
                       int           Q,
                       int           T );

/*
 *  FUNCTION: EDGEBOUNDS_Copy()
 *  SYNOPSIS: Create a deep copy of <edg_src> and store it in <edg_dest>.
 */
EDGEBOUNDS* EDGEBOUNDS_Copy(  EDGEBOUNDS*         edg_dest,
                              const EDGEBOUNDS*   edg_src );

/*
 *  FUNCTION: EDGEBOUNDS_Get()
 *  SYNOPSIS: Return pointer to BOUND at index <i>.
 */
BOUND* EDGEBOUNDS_Get( EDGEBOUNDS*   edg,
                       int           i );

/*
 *  FUNCTION: EDGEBOUNDS_Pushback()
 *  SYNOPSIS: Add BOUND to EDGEBOUNDS list.
 */
void EDGEBOUNDS_Pushback( EDGEBOUNDS*  edg,
                          BOUND*       bnd );

/*
 *  FUNCTION: EDGEBOUNDS_Pushback_Head()
 *  SYNOPSIS: Add head index and row/diag id to lists.
 */
void EDGEBOUNDS_Pushback_Head( EDGEBOUNDS* edg,
                               int         id,
                               int         head );

/*
 *  FUNCTION: EDGEBOUNDS_Insert()
 *  SYNOPSIS: Insert/Overwrite bound into <i> index of Edgebound list.
 */
void EDGEBOUNDS_Insert( EDGEBOUNDS*    edg,
                        int            i,
                        BOUND*         bnd );

/*
 *  FUNCTION: EDGEBOUNDS_Delete()
 *  SYNOPSIS: Delete BOUND at <i> index and fill from end of list.
 */
void EDGEBOUNDS_Delete( EDGEBOUNDS*  edg,
                        int          i );

/*
 *  FUNCTION: EDGEBOUNDS_Clear()
 *  SYNOPSIS: Remove all BOUNDS from EDGEBOUND list.
 */
void EDGEBOUNDS_Clear( EDGEBOUNDS* edg );

/*
 *  FUNCTION: EDGEBOUNDS_Resize()
 *  SYNOPSIS: Resize number of BOUNDS allocated in EDGEBOUND object (does not downsize).
 */
void EDGEBOUNDS_GrowTo( EDGEBOUNDS* edg,
                        int         size );

/*
 *  FUNCTION: EDGEBOUNDS_Resize()
 *  SYNOPSIS: Resize number of BOUNDS in EDGEBOUND object.
 */
void EDGEBOUNDS_Resize(EDGEBOUNDS* edg,
                       int         size);

/*
 *  FUNCTION:  EDGEBOUNDS_Reverse()
 *  SYNOPSIS:  Reverse order of edgebound list.
 */
void EDGEBOUNDS_Reverse(EDGEBOUNDS *edg);

/*
 *  FUNCTION:  EDGEBOUNDS_Index()
 *  SYNOPSIS:  Index locations in EDGEBOUND list that start each unique BOUND id.
 *             Assumes <edg> is sorted.
 */
void EDGEBOUNDS_Index(EDGEBOUNDS *edg);

/*
 *  FUNCTION:  EDGEBOUNDS_Sort()
 *  SYNOPSIS:  Sort <edg> bound list by id, lb, rb.
 */
void EDGEBOUNDS_Sort( EDGEBOUNDS*   edg );

/*
 *  FUNCTION:  EDGEBOUNDS_Merge()
 *  SYNOPSIS:  Merge <edg> bound list by combining overlapping ranges.
 *             Assumes that <edg> is sorted.
 */
void EDGEBOUNDS_Merge( EDGEBOUNDS*   edg );

/*
 *  FUNCTION:  EDGEBOUNDS_Count_Cells()
 *  SYNOPSIS:  Count the number of cells covered by <edg>.
 *             Assumes <edg> is sorted and merged.
 */
int EDGEBOUNDS_Count_Cells( EDGEBOUNDS*   edg );

/*
 *  FUNCTION: EDGEBOUNDS_Print()
 *  SYNOPSIS: Print EDGEBOUND object to file.
 */
void EDGEBOUNDS_Dump(EDGEBOUNDS* edg,
                     FILE*       fp);

/*
 *  FUNCTION: EDGEBOUNDS_Dump()
 *  SYNOPSIS: Output EDGEBOUND object to file.
 */
void EDGEBOUNDS_Save(EDGEBOUNDS*  edg,
                     const char*  _filename_);

/*
 *  FUNCTION: EDGEBOUNDS_Compare()
 *  SYNOPSIS: Compare two EDGEBOUNDS objects.  Return 0 if equal.
 */
int EDGEBOUNDS_Compare( EDGEBOUNDS*    edg_a,
                        EDGEBOUNDS*    edg_b );

/*
 *  FUNCTION: EDGEBOUNDS_Count()
 *  SYNOPSIS: Count the number of cells in edgebound.
 */
int EDGEBOUNDS_Count(EDGEBOUNDS*    edg);

/*
 *  FUNCTION: EDGEBOUNDS_Validate()
 *  SYNOPSIS: Verifies that edgebound ranges don't go out-of-bounds of containing matrix dimensions.
 */
int EDGEBOUNDS_Validate(EDGEBOUNDS *edg);

#endif /* _EDGEBOUND_H */