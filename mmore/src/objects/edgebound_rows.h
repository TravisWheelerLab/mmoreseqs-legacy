/*******************************************************************************
 *  FILE:      edgebound.c
 *  PURPOSE:   EDGEBOUND_ROWS Object.
 *             For building EDGEBOUNDS on the fly.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _EDGEBOUND_ROWS_H
#define _EDGEBOUND_ROWS_H

/*! FUNCTION:  EDGEBOUND_ROWS_Create()
 *  SYNOPSIS:  Create new EDGEBOUND_ROWS object and returns pointer.
 */
EDGEBOUND_ROWS* 
EDGEBOUND_ROWS_Create();

/*! FUNCTION:  EDGEBOUNDS_Create()
 *  SYNOPSIS:  Create new EDGEBOUNDS object with chosen size and returns pointer.
 */
EDGEBOUND_ROWS* 
EDGEBOUND_ROWS_Create_by_Size(   int   Q, 
                                 int   T );

/*! FUNCTION: EDGEBOUND_ROWS_Destroy()
 *  SYNOPSIS: Frees all memory from EDGEBOUND_ROWS object.
 */
EDGEBOUND_ROWS* 
EDGEBOUND_ROWS_Destroy( EDGEBOUND_ROWS*  edg );

/*! FUNCTION: EDGEBOUND_ROWS_Reuse()
 *  SYNOPSIS: Reuses EDGEBOUND_ROWS by "clearing" edgebound list (does not downsize).
 */
void 
EDGEBOUND_ROWS_Reuse(   EDGEBOUND_ROWS*   edg,
                        int               Q,
                        int               T );

/*! FUNCTION: EDGEBOUND_ROWS_Clear()
 *  SYNOPSIS: Reuses EDGEBOUND_ROWS by "clearing" edgebound list (does not downsize).
 */
void 
EDGEBOUND_ROWS_Clear( EDGEBOUND_ROWS*   edg );

/*! FUNCTION: EDGEBOUND_ROWS_GrowTo()
 *  SYNOPSIS: Resizes EDGEBOUND_ROWS if new size exceeds old size
 */
void 
EDGEBOUND_ROWS_GrowTo(  EDGEBOUND_ROWS*  edg,
                        int              size );

/*! FUNCTION: EDGEBOUND_ROWS_Resize()
 *  SYNOPSIS: Resizes EDGEBOUND_ROWS if new size exceeds old size
 */
void 
EDGEBOUND_ROWS_Resize(  EDGEBOUND_ROWS*  edg,
                        int              size );

/*! FUNCTION: EDGEBOUND_ROWS_Get_RowSize()
 *  SYNOPSIS: Get the size of row <q_0>.
 */
int 
EDGEBOUND_ROWS_Get_RowSize(   EDGEBOUND_ROWS*   edg,
                              int               q_0 );

/*! FUNCTION: EDGEBOUND_ROWS_Get()
 *  SYNOPSIS: Return pointer to EDGEBOUND for absolute index <i>.
 */
BOUND* 
EDGEBOUND_ROWS_Get(  EDGEBOUND_ROWS*   edg,
                     int               i );

/*! FUNCTION: EDGEBOUND_ROWS_Get_byRow()
 *  SYNOPSIS: Get pointer to <i_0>th bound on <q_0>th row.
 */
BOUND* 
EDGEBOUND_ROWS_Get_byRow(  EDGEBOUND_ROWS*      edg,
                           int                  q_0,
                           int                  i_0 );

/*! FUNCTION: EDGEBOUND_ROWS_GetLast_byRow()
 *  SYNOPSIS: Gets pointer to the last bound on <q_0>th row.
 */
BOUND* 
EDGEBOUND_ROWS_GetLast_byRow(    EDGEBOUND_ROWS*      edg,
                                 int                  q_0 );

/*! FUNCTION: EDGEBOUND_ROWS_Pushback()
 *  SYNOPSIS: Add BOUND <bnd> to EDGEBOUND_ROWS list at row index <row_id>.
 */
void 
EDGEBOUND_ROWS_Pushback(   EDGEBOUND_ROWS*   edg,
                           int               row_id,
                           BOUND*            bnd );

/*! FUNCTION: EDGEBOUND_ROWS_Copy()
 *  SYNOPSIS: Create a deep copy of <edg_src> and store it in <edg_dest>.
 */
EDGEBOUND_ROWS* 
EDGEBOUND_ROWS_Copy(    EDGEBOUND_ROWS*    edg_dest,
                        EDGEBOUND_ROWS*    edg_src );

/*! FUNCTION: EDGEBOUND_ROWS_Clear()
 *  SYNOPSIS: Remove all BOUNDS from EDGEBOUND list.
 */
void 
EDGEBOUND_ROWS_Clear( EDGEBOUND_ROWS* edg );

/*! FUNCTION: EDGEBOUND_ROWS_IntegrateDiag_Fwd()
 *  SYNOPSIS: Add antidiagonal bound into row-wise bounds, for the Forward Cloud Search.
 *            Looks at each cell individually in the antidiagonal.
 *            If it is right-side adjacent to the current open bound (within a tolerance value), it extends it.
 *            Otherwise, it creates a new edgebound and adds it to the list.
 *            Edgebound lists sizes are determined at compile time and do not resize.
 *            If size is exceeded, program terminates with error.
 */
void 
EDGEBOUND_ROWS_IntegrateDiag_Fwd(   EDGEBOUND_ROWS*   edg,
                                    BOUND*            bnd );

/*! FUNCTION: EDGEBOUND_ROWS_IntegrateDiag_Bck()
 *  SYNOPSIS: Add antidiagonal bound into row-wise bounds, for the Backward Cloud Search.
 *            Looks at each cell individually in the antidiagonal.
 *            If it is left-side adjacent to the current open bound (within a tolerance value), it extends it.
 *            Otherwise, it creates a new edgebound and adds it to the list.
 *            Edgebound lists sizes are determined at compile time and do not resize.
 *            If size is exceeded, program terminates with error.
 */
void 
EDGEBOUND_ROWS_IntegrateDiag_Bck(   EDGEBOUND_ROWS*   edg,
                                    BOUND*            bnd );

/*! FUNCTION: EDGEBOUND_ROWS_Convert()
 *  SYNOPSIS: Convert EDGEBOUNDS_ROWS <edg_in> to EDGEBOUND <edg_out>.
 */
void 
EDGEBOUND_ROWS_Convert(    EDGEBOUND_ROWS*   edg_in,
                           EDGEBOUNDS*       edg_out );

/*! FUNCTION: EDGEBOUND_ROWS_Dump()
 *  SYNOPSIS: Print EDGEBOUND object to file.
 */
void 
EDGEBOUND_ROWS_Dump(    EDGEBOUND_ROWS*   edg,
                        FILE*             fp);

/*! FUNCTION: EDGEBOUND_ROWS_Compare()
 *  SYNOPSIS: Compare two EDGEBOUND_ROWS objects.  Return 0 if equal.
 */
int 
EDGEBOUND_ROWS_Compare(    EDGEBOUND_ROWS*    edg_a,
                           EDGEBOUND_ROWS*    edg_b );

#endif /* _EDGEBOUND_ROWS_H */