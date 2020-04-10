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
void F_INDEX_Destroy(F_INDEX* index);

/* Add INDEX_NODE to F_INDEX nodes array */
void F_INDEX_PushBack(F_INDEX*      index,
                      F_INDEX_NODE  node);

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

/* binary search for node in array in F_INDEX */
int F_INDEX_Search(F_INDEX* index,
                   char*    search_term);

/* sends F_INDEX data to output file */
void F_INDEX_Dump(F_INDEX* index,
                  FILE*    fp);

#endif /* _FILE_INDEX_H_ */

