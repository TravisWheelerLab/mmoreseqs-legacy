/*******************************************************************************
 *  FILE:      edgebound_rows.h
 *  PURPOSE:   EDGEROWS Object.
 *             For building EDGEBOUNDS on the fly.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _EDGEROWS_H
#define _EDGEROWS_H

/*! FUNCTION:  EDGEROWS_Create()
 *  SYNOPSIS:  Create new EDGEROWS object and returns pointer.
 */
EDGEROWS* 
EDGEROWS_Create();

/*! FUNCTION:  EDGEBOUNDS_Create_by_Size()
 *  SYNOPSIS:  Create new EDGEBOUNDS object with chosen size and returns pointer.
 *             Caller must call EDGEROWS_Reuse() before use.
 */
EDGEROWS* 
EDGEROWS_Create_by_Size(   int      Q,
                                 int      T );

/*! FUNCTION: EDGEROWS_Destroy()
 *  SYNOPSIS: Frees all memory from EDGEROWS object.
 */
EDGEROWS* 
EDGEROWS_Destroy( EDGEROWS*  edg );

/*! FUNCTION: EDGEROWS_Reuse()
 *  SYNOPSIS: Reuses EDGEROWS by resizing if too a"clearing" edgebound list (does not downsize).
 */
void 
EDGEROWS_Reuse(   EDGEROWS*   edg,
                        int               Q,
                        int               T,
                        RANGE             Q_range );

/*! FUNCTION: EDGEROWS_Clear()
 *  SYNOPSIS: Reuses EDGEROWS by "clearing" edgebound list (does not downsize).
 */
void 
EDGEROWS_Clear( EDGEROWS*   edg );

/*! FUNCTION: EDGEROWS_GrowTo()
 *  SYNOPSIS: Resizes EDGEROWS if new size exceeds old size
 */
void 
EDGEROWS_GrowTo(  EDGEROWS*  edg,
                        int              size );

/*! FUNCTION: EDGEROWS_Resize()
 *  SYNOPSIS: Resizes EDGEROWS if new size exceeds old size
 */
void 
EDGEROWS_Resize(  EDGEROWS*  edg,
                        int              size );

/*! FUNCTION: EDGEROWS_Get_RowSize()
 *  SYNOPSIS: Get the size of row <q_0>.
 */
int 
EDGEROWS_Get_RowSize(   EDGEROWS*   edg,
                              int               q_0 );

/*! FUNCTION: EDGEROWS_Get()
 *  SYNOPSIS: Return pointer to EDGEBOUND for absolute index <i>.
 */
BOUND* 
EDGEROWS_Get(  EDGEROWS*   edg,
                     int               i );

/*! FUNCTION: EDGEROWS_Get_byRow()
 *  SYNOPSIS: Get pointer to <i_0>th bound on <q_0>th row.
 */
BOUND* 
EDGEROWS_Get_byRow(  EDGEROWS*      edg,
                           int                  q_0,
                           int                  i_0 );

/*! FUNCTION: EDGEROWS_GetLast_byRow()
 *  SYNOPSIS: Gets pointer to the last bound on <q_0>th row.
 */
BOUND* 
EDGEROWS_GetLast_byRow(    EDGEROWS*      edg,
                                 int                  q_0 );

/*! FUNCTION: EDGEROWS_Pushback()
 *  SYNOPSIS: Add BOUND <bnd> to EDGEROWS list at row index <row_id>.
 */
void 
EDGEROWS_Pushback(   EDGEROWS*   edg,
                           int               row_id,
                           BOUND*            bnd );

/*! FUNCTION: EDGEROWS_Copy()
 *  SYNOPSIS: Create a deep copy of <edg_src> and store it in <edg_dest>.
 */
EDGEROWS* 
EDGEROWS_Copy(    EDGEROWS*    edg_dest,
                        EDGEROWS*    edg_src );

/*! FUNCTION: EDGEROWS_Clear()
 *  SYNOPSIS: Remove all BOUNDS from EDGEBOUND list.
 */
void 
EDGEROWS_Clear( EDGEROWS* edg );

/*! FUNCTION: EDGEROWS_IntegrateDiag_Fwd()
 *  SYNOPSIS: Add antidiagonal bound into row-wise bounds, for the Forward Cloud Search.
 *            Looks at each cell individually in the antidiagonal.
 *            If it is right-side adjacent to the current open bound (within a tolerance value), it extends it.
 *            Otherwise, it creates a new edgebound and adds it to the list.
 *            Edgebound lists sizes are determined at compile time and do not resize.
 *            If size is exceeded, program terminates with error.
 */
void 
EDGEROWS_IntegrateDiag_Fwd(   EDGEROWS*   edg,
                                    BOUND*            bnd );

/*! FUNCTION: EDGEROWS_IntegrateDiag_Bck()
 *  SYNOPSIS: Add antidiagonal bound into row-wise bounds, for the Backward Cloud Search.
 *            Looks at each cell individually in the antidiagonal.
 *            If it is left-side adjacent to the current open bound (within a tolerance value), it extends it.
 *            Otherwise, it creates a new edgebound and adds it to the list.
 *            Edgebound lists sizes are determined at compile time and do not resize.
 *            If size is exceeded, program terminates with error.
 */
void 
EDGEROWS_IntegrateDiag_Bck(   EDGEROWS*   edg,
                                    BOUND*            bnd );

/*! FUNCTION: EDGEROWS_Convert()
 *  SYNOPSIS: Convert EDGEBOUNDS_ROWS <edg_in> to EDGEBOUND <edg_out>.
 */
void 
EDGEROWS_Convert(    EDGEROWS*   edg_in,
                           EDGEBOUNDS*       edg_out );

/*! FUNCTION: EDGEROWS_Dump()
 *  SYNOPSIS: Print EDGEBOUND object to file.
 */
void 
EDGEROWS_Dump(    EDGEROWS*   edg,
                        FILE*             fp);

/*! FUNCTION: EDGEROWS_Compare()
 *  SYNOPSIS: Compare two EDGEROWS objects.  Return 0 if equal.
 */
int 
EDGEROWS_Compare(    EDGEROWS*    edg_a,
                           EDGEROWS*    edg_b );

#endif /* _EDGEROWS_H */