/*******************************************************************************
 *  FILE:      f_index.h
 *  PURPOSE:   FILE_INDEX Object.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _FILE_INDEX_H_
#define _FILE_INDEX_H_

/*!  FUNCTION:    F_INDEX_Create()
 *   SYNOPSIS:    Creates an instance of F_INDEX.
 *                Initial node length set to min_size.  
 */
F_INDEX*  
F_INDEX_Create();

/*!  FUNCTION:    F_INDEX_Create()
 *   SYNOPSIS:    Destroys instance of F_INDEX and frees memory.
 */
F_INDEX* 
F_INDEX_Destroy( F_INDEX* index );

/*!  FUNCTION:    F_INDEX_Reuse()
 *   SYNOPSIS:    Reuse an instance of F_INDEX.
 */
void 
F_INDEX_Reuse( F_INDEX* index );

/*!  FUNCTION:    F_INDEX_Pushback()
 *   SYNOPSIS:    Add F_INDEX_NODE to F_INDEX node vector.  
 *                Resizes vector if necessary.
 */
void 
F_INDEX_Pushback(    F_INDEX*       index,
                     F_INDEX_NODE*  node );

/*!  FUNCTION:    F_INDEX_Resize()
 *   SYNOPSIS:    Resizes the length of the F_INDEX nodes array.
 */
void 
F_INDEX_Resize(   F_INDEX*    index,
                  int         size );

/*!  FUNCTION:    F_INDEX_Get()
 *   SYNOPSIS:    Gets the NODE position <idx> in <index>.
 */
F_INDEX_NODE* 
F_INDEX_Get(   F_INDEX*    index,
               int         idx );

/*!  FUNCTION:    F_INDEX_Sort_by_Name()
 *   SYNOPSIS:    Sorts F_INDEX nodes by name.
 */
void 
F_INDEX_Sort_by_Name( F_INDEX*    index );

/*!  FUNCTION:    F_INDEX_Sort_by_Id()
 *   SYNOPSIS:    Sorts F_INDEX nodes by name.
 */
void 
F_INDEX_Sort_by_Id(   F_INDEX*    index );

/*!  FUNCTION:    F_INDEX_Quikort()
 *   SYNOPSIS:    Recursive quicksort of F_INDEX node subarray on range (lo,hi). 
 */
void 
F_INDEX_Quiksort(    F_INDEX_NODE*  arr,    /* F_INDEX node array to be sorted */
                     int            lo,     /* lower end of range in subarray  */
                     int            hi );   /* upper end of range in subarray  */

/*  FUNCTION:    F_INDEX_Swap()
 *  SYNOPSIS:    Swap the ith and jth node in the array.
 */
void 
F_INDEX_Swap(     F_INDEX_NODE*  arr,
                  int            i,
                  int            j );

/*!  FUNCTION:    F_INDEX_Compare_by_Name()
 *   SYNOPSIS:    Compare <a> and <b> node in the array by NAME.
 */
int 
F_INDEX_Compare_by_Name(   const void*    a,
                           const void*    b );

/*!  FUNCTION:    F_INDEX_Compare_by_Id()
 *   SYNOPSIS:    Compare <a> and <b> node in the array by ID.
 */
int 
F_INDEX_Compare_by_Id(  const void*  a,
                        const void*  b );

/*!  FUNCTION:   F_INDEX_Getby_Name()
 *   SYNOPSIS:   Runs NAME Search for <search_term> and returns node.
 *     RETURN:   Search result; NULL if no result found.
 */
F_INDEX_NODE* 
F_INDEX_Getby_Name(    F_INDEX*    index, 
                        char*       search_term );

/*!  FUNCTION:   F_INDEX_Search_Name()
 *   SYNOPSIS:   Binary search (by name) for node in array in F_INDEX.
 *               Assumes F_INDEX is sorted by Name.
 *     RETURN:   index of search result; -1 if no result found.
 */
int 
F_INDEX_Search_Name(    F_INDEX*    index,
                        char*       search_term );

/*!  FUNCTION:   F_INDEX_Getby_Id()
 *   SYNOPSIS:   Runs ID Search for <search_term> and returns node.
 *     RETURN:   Search result; NULL if no result found.
 */
F_INDEX_NODE* 
F_INDEX_Getby_Id(   F_INDEX* index, 
                     int      search_term );

/*!  FUNCTION:    F_INDEX_Search_Id()
 *   SYNOPSIS:    Binary search (by id) for node in array in F_INDEX. 
 *                Assumes F_INDEX is sorted by Id.
 *     RETURN:    index of search result; -1 if no result found.
 */
int 
F_INDEX_Search_Id(   F_INDEX*    index,
                     int         search_term );

/*!  FUNCTION:    F_INDEX_Save()
 *   SYNOPSIS:    Save F_INDEX data to filepath.
 */
void 
F_INDEX_Save(  F_INDEX*   index,
               char*      _filename_ );

/*!  FUNCTION:    F_INDEX_Dump()
 *   SYNOPSIS:    Send F_INDEX data to file.
 */
void 
F_INDEX_Dump(  F_INDEX*   index,
               FILE*      fp );

/*!  FUNCTION:    F_INDEX_Node_Dump()
 *   SYNOPSIS:    Output F_INDEX node data to file pointer at index.
 */
void 
F_INDEX_Node_Dump(   F_INDEX*    index,
                     int         id,
                     FILE*       fp );

/*!  FUNCTION:    F_INDEX_UnitTest
 *   SYNOPSIS:    Unit Test for F_INDEX.
 */
void 
F_INDEX_UnitTest();

#endif /* _FILE_INDEX_H_ */

