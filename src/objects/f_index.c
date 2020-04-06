/*******************************************************************************
 *  FILE:      f_index.h
 *  PURPOSE:   F_INDEX Object.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* local imports */
#include "objects/structs.h"
#include "utilities/utility.h"

/* unit test imports */
#include "index_parser.h"

/* header */
#include "f_index.h"

/* === PRIVATE FUNCTIONS === */
void F_INDEX_Resize(F_INDEX* index,
                    int      size);

void F_INDEX_Quiksort(F_INDEX_NODE*  arr,
                      int            lo,
                      int            hi);


/* *******************************************************************
 *    FUNC:    F_INDEX_Create()
 *    DESC:    Creates an instance of F_INDEX.
 *             Initial node length set to min_size.  
 *
 *    ARGS:       <pathname>     Path to file being indexed
 *
 *    RETURN:     Pointer to new F_INDEX object
 * *******************************************************************/
F_INDEX*  F_INDEX_Create(char* pathname)
{
   F_INDEX*   index    = NULL;
   const int  min_size = 16;

   index = (F_INDEX*) malloc( sizeof(F_INDEX) );
   if (index == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc F_INDEX.\n");
      exit(EXIT_FAILURE);
   }

   index->N          = 0;
   index->Nalloc     = min_size;
   index->filepath   = NULL;
   index->nodes      = NULL;
   index->isSorted   = false;

   index->nodes      = (F_INDEX_NODE*) malloc( sizeof(F_INDEX_NODE) * min_size );
   if (index == NULL) {
      fprintf(stderr, "ERROR: Unable to malloc F_INDEX_NODE array for F_INDEX.\n");
      exit(EXIT_FAILURE);
   }

   index->filepath   = strdup(pathname);

   return index;
}

/* *******************************************************************
 *    FUNC:    F_INDEX_Create()
 *    DESC:    Destroys instance of F_INDEX and frees memory.
 *
 *    ARGS:       <index>     F_INDEX object to be freed (cannot be NULL)
 *
 *    RETURN:     None.
 * *******************************************************************/
void F_INDEX_Destroy(F_INDEX* index)
{
   for (int i = 0; i < index->N; i++) {
      free(index->nodes[i].name);
   }
   free(index->nodes);
   free(index);
}

/* *******************************************************************
 *    FUNC:    F_INDEX_PushBack()
 *    DESC:    Add F_INDEX_NODE to F_INDEX node vector.  
 *             Resizes vector if necessary.
 *
 *    ARGS:       <index>     F_INDEX 
 *                <node>      F_INDEX_NODE to be added to F_INDEX
 *
 *    RETURN:     None.
 * *******************************************************************/
void F_INDEX_PushBack(F_INDEX*      index,
                      F_INDEX_NODE  node)
{
   index->nodes[index->N] = node;
   /* allocate space for string */
   index->nodes[index->N].name = strdup(node.name);
   
   index->N++;
   if (index->N >= index->Nalloc) {
      F_INDEX_Resize(index, index->N * 2);
   }
}

/* *******************************************************************
 *    FUNC:    F_INDEX_Resize()
 *    DESC:    Resizes the length of the F_INDEX nodes array.
 *
 *    ARGS:       <index>     F_INDEX 
 *                <size>      Length F_INDEX node array is to be resized to.
 *
 *    RETURN:     None.
 * *******************************************************************/
void F_INDEX_Resize(F_INDEX* index,
                    int      size)
{
   index->Nalloc = size;
   index->nodes  = (F_INDEX_NODE*) realloc(index->nodes, sizeof(F_INDEX_NODE) * size );
   if ( index->nodes == NULL ) {
      fprintf(stderr, "ERROR: Unable to realloc NODES for F_INDEX.\n");
      exit(EXIT_FAILURE);
   }
}

/* *******************************************************************
 *    FUNC:    F_INDEX_Sort()
 *    DESC:    Sorts F_INDEX nodes by name.
 *
 *    ARGS:       <index>     F_INDEX
 *
 *    RETURN:     None.
 * *******************************************************************/
void F_INDEX_Sort(F_INDEX* index)
{
   F_INDEX_Quiksort(index->nodes, 0, index->N);
}

/* *******************************************************************
 *    FUNC:    F_INDEX_Quikort()
 *    DESC:    Recursive quicksort of F_INDEX node subarray on range (lo,hi). 
 *
 *    ARGS:       <arr>       F_INDEX node array to be sorted
 *                <lo>        lower end of range in subarray 
 *                <hi>        upper end of range in subarray          
 *
 *    RETURN:     None.
 * *******************************************************************/
void F_INDEX_Quiksort(F_INDEX_NODE*  arr,
                      int            lo,
                      int            hi)
{
   if (hi - lo <= 1) return;

   /* set pivot element and  */
   int pivot;
   srand(time(NULL));
   pivot = ( rand() % (hi - lo) ) + lo;
   F_INDEX_Swap(arr, 0, pivot);

   int left  = lo;
   int right = hi - 1; 

   /* partition array according to pivot element */
   while (true) {
      /* shift left pointer to the right until a node is found greater than pivot */
      while ( left <= right ) {
         if ( F_INDEX_Compare( arr, left, 0 ) > 0 ) break;
         left++;
      }

      /* shift right pointer to the left until a node is found greater than pivot */
      while ( left <= right ) {
         if ( F_INDEX_Compare( arr, right, 0 ) <= 0 ) break;
         right--;
      }

      /* if left and right indexes have crossed, then partition is done */
      if (left >= right) break;

      F_INDEX_Swap( arr, left, right );   
   }
   /* finally, swap pivot to right-edge of the left partition */
   F_INDEX_Swap( arr, 0, left - 1 );

   /* recurse on subproblems */
   F_INDEX_Quiksort( arr, 0, left - 1 );
   F_INDEX_Quiksort( arr, left, hi );
}

/* swap ith and jth nodes in nodes array */
inline
void F_INDEX_Swap(F_INDEX_NODE*  arr,
                  int            i,
                  int            j)
{
   F_INDEX_NODE tmp;
   tmp    = arr[i];
   arr[i] = arr[j];
   arr[j] = tmp;
}

/* compare ith to jth node in nodes array */
inline
int F_INDEX_Compare(F_INDEX_NODE*  arr,
                    int            i,
                    int            j)
{
   // printf("COMPARE: {%s} VS {%s}\n", arr[i].name, arr[j].name);
   int cmp = strcmp(arr[i].name, arr[j].name);
   return cmp;
}

/* binary search for node in array in F_INDEX */
int F_INDEX_Search(F_INDEX* index,
                   char*    search_term)
{
   int lo  = 0;
   int mid = 0;
   int hi  = index->N;
   int cmp = 0;

   printf("SEARCH TERM: '%s'\n", search_term);
   while (lo <= hi)
   {
      mid = (lo+hi)/2;
      cmp = strcmp( search_term, index->nodes[mid].name );

      if ( cmp < 0 ) {
         hi = mid - 1;
      }
      else 
      if ( cmp > 0 ) {
         lo = mid + 1;
      }
      else {
         return mid;
      }
      // printf("SEARCHING, range=(%d,%d), mid=%d, cmp=%d, term='%s'...\n", lo, hi, mid, cmp, index->nodes[mid].name);
   }
   return -1;
}

/* sends F_INDEX data to output file */
void F_INDEX_Dump(F_INDEX* index,
                  FILE*    fp)
{
   if (fp == NULL) {
      fprintf(stderr, "ERROR: Unable to open file.\n");
      exit(EXIT_FAILURE);
   }

   fprintf(fp, "# === FILE INDEX === #\n");
   fprintf(fp, "%s\t%s\n", "FILE_NAME:",     index->filepath);
   fprintf(fp, "%s\t%d\n", "NUMBER_SEQS:",   index->N);
   fprintf(fp, "# {ID}\t{OFF}\t{NAME}\n");

   for (int i = 0; i < index->N; i++)
   {
      F_INDEX_NODE node = index->nodes[i];
      fprintf(fp, "%d\t%ld\t%s\n", i, node.offset, node.name);
   }
   fprintf(fp, "\n");
}

/* unit test for F_INDEX */
void F_INDEX_UnitTest()
{
   F_INDEX* f_index = F_INDEX_Fasta_Build("src/test.txt");
   F_INDEX_Sort(f_index);
   F_INDEX_Dump(f_index, stdout);

   for (int i = 0; i < f_index->N; i++) 
   {
      char* search_term = f_index->nodes[i].name;
      int   find        = F_INDEX_Search(f_index, search_term);
      printf("INDEX: %d, FIND: %s\n", find, search_term);
   }
   exit(EXIT_SUCCESS);
}