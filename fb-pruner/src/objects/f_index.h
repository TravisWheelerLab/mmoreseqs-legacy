/*******************************************************************************
 *  FILE:      f_index.h
 *  PURPOSE:   FILE_INDEX Object.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

#ifndef _FILE_INDEX_H_
#define _FILE_INDEX_H_

/* Constructor */
F_INDEX*  F_INDEX_Create();

/* Destructor */
void* F_INDEX_Destroy(F_INDEX* index);

/*
 *    FUNC:    F_INDEX_Reuse()
 *    DESC:    Reuse an instance of F_INDEX.
 */
void F_INDEX_Reuse( F_INDEX* index );

/* Add INDEX_NODE to F_INDEX nodes array */
void F_INDEX_Pushback(F_INDEX*      index,
                      F_INDEX_NODE* node);

/* Resizes the length of the nodes array */
void F_INDEX_Resize(F_INDEX* index,
                    int      size);

/* sorts index by name */
void F_INDEX_Sort_by_Name(F_INDEX* index);

/* sorts index by name */
void F_INDEX_Sort_by_Id(F_INDEX* index);

/* run quiksort on node array of given length */
void F_INDEX_Quiksort(F_INDEX_NODE*  arr,
                      int            lo,
                      int            hi);

/* swap ith and jth nodes in nodes array */
void F_INDEX_Swap(F_INDEX_NODE*  arr,
                  int            i,
                  int            j);

/* compare ith to jth node in nodes array */
int F_INDEX_Compare_by_Name(const void* a,
                            const void* b);

/* compare ith to jth node in nodes array */
int F_INDEX_Compare_by_Id(const void*  a,
                          const void*  b);

/* *******************************************************************
 * FUNCTION:   F_INDEX_Search_Name()
 * SYNOPSIS:   Binary search (by name) for node in array in F_INDEX.
 *             Assumes F_INDEX is sorted by Name.
 * RETURN:     index of search result; -1 if no result found.
 * *******************************************************************/
int F_INDEX_Search_Name(F_INDEX* index,
                        char*    search_term);

/* *******************************************************************
 *    FUNC:    F_INDEX_Search_Id()
 *    DESC:    Binary search (by id) for node in array in F_INDEX. 
 *             Assumes F_INDEX is sorted by Id.
 * RETURN:     index of search result; -1 if no result found.
 * *******************************************************************/
int F_INDEX_Search_Id( F_INDEX* index,
                       int      search_term);

/*
 *    FUNCTION:    F_INDEX_Save()
 *    SYNOPSIS:    Save F_INDEX data to filepath.
 */
void F_INDEX_Save( F_INDEX*   index,
                   char*      _filename_ );

/* sends F_INDEX data to output file */
void F_INDEX_Dump(F_INDEX* index,
                  FILE*    fp);

/* *******************************************************************
 *    FUNC:    F_INDEX_Node_Dump()
 *    DESC:    Output F_INDEX node data to file pointer at index.
 * *******************************************************************/
void F_INDEX_Node_Dump( F_INDEX*    index,
                        int         id,
                        FILE*       fp );

#endif /* _FILE_INDEX_H_ */

