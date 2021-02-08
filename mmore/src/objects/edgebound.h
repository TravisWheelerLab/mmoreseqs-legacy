/*******************************************************************************
 *  FILE:      edgebound.c
 *  PURPOSE:   EDGEBOUNDS Object
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _EDGEBOUND_H
#define _EDGEBOUND_H

/*! FUNCTION:  EDGEBOUNDS_Create()
 *  SYNOPSIS:  Create new EDGEBOUNDS object and returns pointer.
 */
EDGEBOUNDS* 
EDGEBOUNDS_Create();

/*! FUNCTION:  EDGEBOUNDS_Create_by_Size()
 *  SYNOPSIS:  Create new EDGEBOUNDS object with chosen size and returns pointer.
 */
EDGEBOUNDS* 
EDGEBOUNDS_Create_by_Size( const int size );

/*! FUNCTION: EDGEBOUNDS_Destroy()
 *  SYNOPSIS: Frees all memory from EDGEBOUNDS object.
 */
EDGEBOUNDS* 
EDGEBOUNDS_Destroy( EDGEBOUNDS*  edg );

/*! FUNCTION: EDGEBOUNDS_Reuse()
 *  SYNOPSIS: Reuses EDGEBOUNDS by "clearing" edgebound list (does not realloc).
 */
void EDGEBOUNDS_Reuse( EDGEBOUNDS*   edg, 
                       int           Q,
                       int           T );

/*! FUNCTION: EDGEBOUNDS_Copy()
 *  SYNOPSIS: Create a deep copy of <edg_src> and store it in <edg_dest>.
 */
EDGEBOUNDS* 
EDGEBOUNDS_Copy(    EDGEBOUNDS*          edg_dest,
                    const EDGEBOUNDS*    edg_src );

/*! FUNCTION: EDGEBOUNDS_Get()
 *  SYNOPSIS: Return pointer to BOUND at index <i>.
 */
BOUND* 
EDGEBOUNDS_Get(   EDGEBOUNDS*   edg,
                  int           i );

/*! FUNCTION: EDGEBOUNDS_GetSize()
 *  SYNOPSIS: Get length of <edg>.
 */
int 
EDGEBOUNDS_GetSize(  EDGEBOUNDS*   edg );

/*! FUNCTION: EDGEBOUNDS_SetSize()
 *  SYNOPSIS: Set length of <edg> to <size>.
 */
STATUS_FLAG 
EDGEBOUNDS_SetSize(    EDGEBOUNDS*    edg,
                        int            size );

/*! FUNCTION:  EDGEBOUNDS_Search()
 *  SYNOPSIS:  Binary search edgebounds for bound containing cell (q_0, t_0).
 *             Assumes edgebounds are sorted and merged.
 *  RETURN:    Return index of edgebound, or -1 if not contained.
 */
int 
EDGEBOUNDS_Search(   EDGEBOUNDS*    edg,     /* edgebounds  */
                     int            q_0,     /* row/diag index, position in query */
                     int            t_0 );   /* column index, position in target */

/*! FUNCTION: EDGEBOUNDS_Pushback()
 *  SYNOPSIS: Add BOUND to EDGEBOUNDS list.
 */
STATUS_FLAG 
EDGEBOUNDS_Pushback(    EDGEBOUNDS*  edg,
                        BOUND*       bnd );

/*! FUNCTION: EDGEBOUNDS_Insert()
 *  SYNOPSIS: Insert/Overwrite bound into <i> index of Edgebound list.
 */
STATUS_FLAG 
EDGEBOUNDS_Insert(   EDGEBOUNDS*    edg,
                     int            i,
                     BOUND*         bnd );

/*! FUNCTION: EDGEBOUNDS_Delete()
 *  SYNOPSIS: Delete BOUND at <i> index and fill from end of list <N-1>, then decrement list size.
 */
STATUS_FLAG 
EDGEBOUNDS_Delete(   EDGEBOUNDS*    edg,
                     int            i );

/*! FUNCTION: EDGEBOUNDS_Resize()
 *  SYNOPSIS: Resize number of BOUNDS allocated in EDGEBOUND object (does not downsize).
 */
STATUS_FLAG 
EDGEBOUNDS_GrowTo(   EDGEBOUNDS*    edg,
                     int            size );

/*! FUNCTION: EDGEBOUNDS_Clear()
 *  SYNOPSIS: Remove all BOUNDS from EDGEBOUND list.
 */
void 
EDGEBOUNDS_Clear( EDGEBOUNDS*  edg );

/*! FUNCTION: EDGEBOUNDS_Resize()
 *  SYNOPSIS: Resize number of BOUNDS in EDGEBOUND object.
 */
void 
EDGEBOUNDS_Resize(   EDGEBOUNDS*     edg,
                     int             size );

/*! FUNCTION:  EDGEBOUNDS_Reverse()
 *  SYNOPSIS:  Reverse order of edgebound list.
 */
STATUS_FLAG 
EDGEBOUNDS_Reverse( EDGEBOUNDS*    edg );

/*! FUNCTION:  EDGEBOUNDS_Index()
 *  SYNOPSIS:  Index locations in EDGEBOUND list that start each unique BOUND id.
 *             Assumes <edg> is sorted.
 */
STATUS_FLAG 
EDGEBOUNDS_Index( EDGEBOUNDS*  edg );

/*! FUNCTION:  EDGEBOUNDS_NxtRow()
 *  SYNOPSIS:  Using iterators, gets the row index range <r_0> that are on <q_0> position.
 *             Skips over rows less than <q_0>.  Presumes that edgebounds is sorted.  
 */
STATUS_FLAG 
EDGEBOUNDS_NxtRow(   EDGEBOUNDS*          edg,     /* edgebounds */
                     int*                 r_0b,    /* row range begin */
                     int*                 r_0e,    /* row range end */
                     int                  q_0 );   /* query sequence position */

/*! FUNCTION:  EDGEBOUNDS_PrvRow()
 *  SYNOPSIS:  Using iterators, gets the row index range <r_0> that are on <q_0> position.
 *             Skips over rows greater than <q_0>.  Presumes that edgebounds is sorted.  
 */
STATUS_FLAG 
EDGEBOUNDS_PrvRow(   EDGEBOUNDS*          edg,     /* edgebounds */
                     int*                 r_0b,    /* row range begin */
                     int*                 r_0e,    /* row range end */
                     int                  q_0 );   /* query sequence position */

/*! FUNCTION:  EDGEBOUNDS_Sort()
 *  SYNOPSIS:  Sort <edg> bound list by id, lb, rb.
 */
void 
EDGEBOUNDS_Sort( EDGEBOUNDS*   edg );

/*! FUNCTION:  EDGEBOUNDS_Sort_Sub()
 *  SYNOPSIS:  Sort the edgebounds subarray on range (beg, end]. Sorts in place.
 */
void 
EDGEBOUNDS_Sort_Sub(    EDGEBOUNDS*    edg,
                        int            beg,
                        int            end );

/*! FUNCTION:  EDGEBOUNDS_Sort_Sub_Selectsort()
 *  SYNOPSIS:  Selection Sorts subarray of <vec> data in ascending order on range (beg,end].  
 */
void 
EDGEBOUNDS_Sort_Sub_Selectsort(     EDGEBOUNDS*    edg,
                                    int            beg,
                                    int            end );

/*! FUNCTION:  EDGEBOUNDS_Sort_Sub_Quicksort()
 *  SYNOPSIS:  Quick Sorts subarray of <vec> data in ascending order on range (beg,end].  
 */
void 
EDGEBOUNDS_Sort_Sub_Quicksort(   EDGEBOUNDS*    edg,
                                 int            beg,
                                 int            end );

/*! FUNCTION:  EDGEBOUNDS_Swap()
 *  SYNOPSIS:  Swaps the values of <vec> at indexes <i> and <j>.
 *             Warning: Only checks for Out-of-Bounds when in DEBUG.
 */
void 
EDGEBOUNDS_Swap(  EDGEBOUNDS*    edg,
                  int            i,
                  int            j );

/*! FUNCTION:  EDGEBOUNDS_Merge()
 *  SYNOPSIS:  Merge <edg> bound list by combining overlapping ranges.
 *             Assumes that <edg> is sorted.
 */
void 
EDGEBOUNDS_Merge( EDGEBOUNDS*   edg );

/*! FUNCTION:  EDGEBOUNDS_Merge_Sub()
 *  SYNOPSIS:  Merge <edg>'s subarray of bound list by combining overlapping ranges. In-place.
 *             Assumes that <edg> is sorted.
 */
void 
EDGEBOUNDS_Merge_Sub(   EDGEBOUNDS*    edg,
                        int            beg,
                        int            end );

/*! FUNCTION: EDGEBOUNDS_Print()
 *  SYNOPSIS: Print EDGEBOUND object to file.
 */
void 
EDGEBOUNDS_Dump(  EDGEBOUNDS*     edg,
                  FILE*           fp );

/*! FUNCTION: EDGEBOUNDS_Print()
 *  SYNOPSIS: Print EDGEBOUND object to file.
 */
void 
EDGEBOUNDS_Sub_Dump(    EDGEBOUNDS*    edg,
                        FILE*          fp,
                        int            beg, 
                        int            end );

/*! FUNCTION: EDGEBOUNDS_Dump()
 *  SYNOPSIS: Output EDGEBOUND object to file.
 */
void 
EDGEBOUNDS_Save(  EDGEBOUNDS*   edg,
                  const char*   _filename_ );

/*! FUNCTION: EDGEBOUNDS_Compare()
 *  SYNOPSIS: Compare two EDGEBOUNDS objects.  Return 0 if equal.
 */
int 
EDGEBOUNDS_Compare(  EDGEBOUNDS*    edg_a,
                     EDGEBOUNDS*    edg_b );

/*! FUNCTION: EDGEBOUNDS_Count()
 *  SYNOPSIS: Count the number of cells in edgebound.
 */
int 
EDGEBOUNDS_Count( EDGEBOUNDS*    edg );

/*! FUNCTION: EDGEBOUNDS_Validate()
 *  SYNOPSIS: Verifies that edgebound ranges don't go out-of-bounds of containing matrix dimensions.
 */
int EDGEBOUNDS_Validate( EDGEBOUNDS*    edg );

/*! FUNCTION: EDGEBOUNDS_Cover_Matrix()
 *  SYNOPSIS: For testing. Creates an edgebounds that covers every cell in DP Matrix with dimensions {Q x T}.
 */
void
EDGEBOUNDS_Cover_Matrix(   EDGEBOUNDS*    edg, 
                           int            Q,
                           int            T );

/*! FUNCTION:  EDGEBOUNDS_Cover_Range()
 *  SYNOPSIS:  For testing.
 *             Creates edgebound space that fills square with Q_range in Query and T_range in Target in DP Matrix.
 */
STATUS_FLAG 
EDGEBOUNDS_Cover_Range(    EDGEBOUNDS*    edg,
                           RANGE          Q_range,
                           RANGE          T_range );

/*! FUNCTION: EDGEBOUNDS_Find_BoundingBox()
 *  SYNOPSIS: Find the min/max range of values contained in the edgebounds.
 *            Assumes edgebounds have been sorted.
 */
int 
EDGEBOUNDS_Find_BoundingBox(  EDGEBOUNDS*   edg,
                              RANGE*        Q_range,
                              RANGE*        T_range );

/*! FUNCTION: EDGEBOUNDS_SetDomain()
 *  SYNOPSIS: Build an EDGEBOUND <edg_out> from QxT EDGEBOUNDS <edg_in>
 *            and constraining the query range to the domain <dom_range> => <q_beg, q_end>.
 *            Simply eliminates query rows outside the range and shifts all query id's by <q_beg>.
 */
int 
EDGEBOUNDS_SetDomain(  EDGEBOUNDS*     edg_in,
                        EDGEBOUNDS*     edg_out,
                        RANGE           dom_range );           

#endif /* _EDGEBOUND_H */