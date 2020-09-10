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
#include "structs.h"
#include "utilities.h"
#include "objects.h"

/* unit test imports */
#include "index_parser.h"

/* header */
#include "f_index.h"

/* === PRIVATE FUNCTIONS === */
void F_INDEX_Resize( F_INDEX*    index,
                     int         size );

void F_INDEX_Quiksort(  F_INDEX_NODE*  arr,
                        int            lo,
                        int            hi );


/*    FUNC:    F_INDEX_Create()
 *    DESC:    Creates an instance of F_INDEX.
 *             Initial node length set to min_size.  
 */
F_INDEX*  F_INDEX_Create()
{
   F_INDEX*   index    = NULL;
   const int  min_size = 32;

   index = (F_INDEX*) ERROR_malloc( sizeof(F_INDEX) );

   index->N             = 0;
   index->Nalloc        = min_size;
   index->nodes         = NULL;

   index->index_path    = NULL;
   index->lookup_path   = NULL;
   index->source_path   = NULL;
   index->delim         = NULL;

   index->sort_type     = SORT_NONE;   /* not sorted */
   index->mmseqs_names  = false;

   index->nodes         = (F_INDEX_NODE*) ERROR_malloc( sizeof(F_INDEX_NODE) * min_size );
   return index;
}

/*    FUNC:    F_INDEX_Create()
 *    DESC:    Destroys instance of F_INDEX and frees memory.
 */
void* F_INDEX_Destroy( F_INDEX* index )
{
   if (index == NULL) return index;

   for (int i = 0; i < index->N; i++) {
      ERROR_free(index->nodes[i].name);
   }
   ERROR_free(index->nodes);

   ERROR_free(index->index_path);
   ERROR_free(index->lookup_path);
   ERROR_free(index->source_path);
   ERROR_free(index->delim);

   ERROR_free(index);
   index = NULL;
   return index;
}

/*    FUNC:    F_INDEX_Reuse()
 *    DESC:    Reuse an instance of F_INDEX.
 */
void F_INDEX_Reuse( F_INDEX* index )
{
   index->N             = 0;

   ERROR_free(index->index_path);
   index->lookup_path   = NULL;
   ERROR_free(index->source_path);
   index->source_path   = NULL;
   ERROR_free(index->delim);
   index->delim         = NULL;

   index->sort_type     = SORT_NONE;
   index->mmseqs_names  = false;
}


/*    FUNC:    F_INDEX_Pushback()
 *    DESC:    Add F_INDEX_NODE to F_INDEX node vector.  
 *             Resizes vector if necessary.
 */
void F_INDEX_Pushback(  F_INDEX*       index,
                        F_INDEX_NODE*  node )
{
   index->nodes[index->N] = *node;
   /* allocate space for string */
   index->nodes[index->N].name = strdup(node->name);
   
   index->N++;
   if (index->N >= index->Nalloc) {
      F_INDEX_Resize(index, index->N * 2);
   }
}

/*    FUNC:    F_INDEX_Resize()
 *    DESC:    Resizes the length of the F_INDEX nodes array.
 */
void F_INDEX_Resize( F_INDEX*    index,
                     int         size )
{
   index->Nalloc = size;
   index->nodes  = (F_INDEX_NODE*) realloc(index->nodes, sizeof(F_INDEX_NODE) * size );
   if ( index->nodes == NULL ) {
      fprintf(stderr, "ERROR: Unable to realloc NODES for F_INDEX.\n");
      exit(EXIT_FAILURE);
   }
}

/* *******************************************************************
 *    FUNC:    F_INDEX_Sort_by_Name()
 *    DESC:    Sorts F_INDEX nodes by name.
 * *******************************************************************/
void F_INDEX_Sort_by_Name( F_INDEX*    index )
{
   // F_INDEX_Quiksort(index->nodes, 0, index->N);
   qsort(index->nodes, index->N, sizeof(F_INDEX_NODE), F_INDEX_Compare_by_Name);
   index->sort_type = SORT_NAME;
}

/* *******************************************************************
 *    FUNC:    F_INDEX_Sort_by_Id()
 *    DESC:    Sorts F_INDEX nodes by name.
 * *******************************************************************/
void F_INDEX_Sort_by_Id(   F_INDEX*    index )
{
   // F_INDEX_Quiksort(index->nodes, 0, index->N);
   qsort(index->nodes, index->N, sizeof(F_INDEX_NODE), F_INDEX_Compare_by_Id);
   index->sort_type = SORT_ID;
}

/* *******************************************************************
 *    FUNC:    F_INDEX_Quikort()
 *    DESC:    Recursive quicksort of F_INDEX node subarray on range (lo,hi). 
 * *******************************************************************/
void F_INDEX_Quiksort(  F_INDEX_NODE*  arr,    /* F_INDEX node array to be sorted */
                        int            lo,     /* lower end of range in subarray  */
                        int            hi )    /* upper end of range in subarray  */
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
         if ( F_INDEX_Compare_by_Name( &(arr[left]), &(arr[0]) ) > 0 ) break;
         left++;
      }

      /* shift right pointer to the left until a node is found greater than pivot */
      while ( left <= right ) {
         if ( F_INDEX_Compare_by_Name( &(arr[right]), &(arr[0]) ) <= 0 ) break;
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

/* *******************************************************************
 *    FUNC:    F_INDEX_Swap()
 *    DESC:    Swap the ith and jth node in the array.
 * *******************************************************************/
inline
void F_INDEX_Swap(   F_INDEX_NODE*  arr,
                     int            i,
                     int            j )
{
   F_INDEX_NODE tmp;
   tmp    = arr[i];
   arr[i] = arr[j];
   arr[j] = tmp;
}

/* *******************************************************************
 *    FUNC:    F_INDEX_Compare_by_Name()
 *    DESC:    Compare <a> and <b> node in the array by NAME.
 * *******************************************************************/
inline
int F_INDEX_Compare_by_Name(  const void* a,
                              const void* b )
{
   F_INDEX_NODE* node_a = (F_INDEX_NODE*)a;
   F_INDEX_NODE* node_b = (F_INDEX_NODE*)b; 
   int cmp = strcmp( node_a->name, node_b->name);
   return cmp;
}

/* *******************************************************************
 *    FUNC:    F_INDEX_Compare_by_Id()
 *    DESC:    Compare <a> and <b> node in the array by ID.
 * *******************************************************************/
inline
int F_INDEX_Compare_by_Id( const void*  a,
                           const void*  b )
{
   F_INDEX_NODE* node_a = (F_INDEX_NODE*)a;
   F_INDEX_NODE* node_b = (F_INDEX_NODE*)b; 
   int cmp = node_a->id - node_b->id;
   return cmp;
}

/*  FUNCTION:   F_INDEX_Search_Name()
 *  SYNOPSIS:   Binary search (by name) for node in array in F_INDEX.
 *              Assumes F_INDEX is sorted by Name.
 *  RETURN:     index of search result; -1 if no result found.
 */
int F_INDEX_Search_Name(   F_INDEX* index,
                           char*    search_term )
{
   int lo  = 0;
   int mid = 0;
   int hi  = index->N;
   int cmp = 0;
   F_INDEX_NODE node;

   #if DEBUG
   {
      if ( index->sort_type != SORT_NAME )
      {
         printf("ERROR: Binary Search of F_INDEX by Name while not sorted by Name.\n");
         exit(EXIT_FAILURE);
      }
   }
   #endif

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
   }
   return -1;
}

F_INDEX_NODE* F_INDEX_Get_by_Id( F_INDEX* index, 
                                 int      search_term )
{
   int res = F_INDEX_Search_Id( index, search_term );
   if (res == -1) return NULL;
   return &(index->nodes[res]);
}

/* *******************************************************************
 *    FUNC:    F_INDEX_Search_Id()
 *    DESC:    Binary search (by id) for node in array in F_INDEX. 
 *             Assumes F_INDEX is sorted by Id.
 * RETURN:     index of search result; -1 if no result found.
 * *******************************************************************/
int F_INDEX_Search_Id( F_INDEX* index,
                       int      search_term )
{
   int lo  = 0;
   int mid = 0;
   int hi  = index->N;
   int cmp = 0;

   #if DEBUG
   {
      if ( index->sort_type != SORT_ID )
      {
         printf("ERROR: Binary Search of F_INDEX by ID while not sorted by ID.\n");
         exit(EXIT_FAILURE);
      }
   }
   #endif
   
   while (lo <= hi)
   {
      mid = (lo+hi)/2;
      cmp = search_term - index->nodes[mid].id;

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
   }
   return -1;
}

/* *******************************************************************
 *    FUNC:    F_INDEX_Save()
 *    DESC:    Save F_INDEX data to filepath.
 * *******************************************************************/
void F_INDEX_Save( F_INDEX*   index,
                   char*      _filename_ )
{
   FILE* fp = fopen( _filename_, "w" );
   F_INDEX_Dump( index, fp );
   printf("F_INDEX saved to: '%s'\n", _filename_);
}

/* *******************************************************************
 *    FUNC:    F_INDEX_Dump()
 *    DESC:    Send F_INDEX data to file.
 * *******************************************************************/
void F_INDEX_Dump( F_INDEX*   index,
                   FILE*      fp )
{
   F_INDEX_NODE*  node;
   char*          name;
   char*          delim = " \t";    /* whitespace */

   if (fp == NULL) {
      fprintf(stderr, "ERROR: Unable to open file.\n" );
      exit(EXIT_FAILURE);
   }

   /* print header (TODO: make option?) */
   int header = false;
   if (header == true) {
      fprintf(fp, "# === FILE INDEX === #\n" );
      // fprintf(fp, "%s\t%s\n", "INDEX_PATH:",    index->index_path );
      fprintf(fp, "%s\t%s\n", "SOURCE_PATH:",   index->source_path );
      fprintf(fp, "%s\t%s\n", "LOOKUP_PATH:",   index->lookup_path );
      fprintf(fp, "%s\t%d\n", "NUMBER_SEQS:",   index->N );
      fprintf(fp, "%s\t%s\n", "MMSEQS_NAMES:",  index->mmseqs_names ? "true" : "false" );
      fprintf(fp, ">{ID}\t{OFF}\t{NAME}\n");
   }
   
   /* print index */
   for (int i = 0; i < index->N; i++)
   {
      node = &(index->nodes[i]);
      name = strtok(node->name, delim);
      fprintf(fp, "%d\t%ld\t%s\t", i, node->offset, name );
      fprintf(fp, "\n");
   }
}

/* *******************************************************************
 *    FUNC:    F_INDEX_Node_Dump()
 *    DESC:    Output F_INDEX node data to file pointer at index.
 * *******************************************************************/
void F_INDEX_Node_Dump( F_INDEX*    index,
                        int         id,
                        FILE*       fp )
{
   F_INDEX_NODE*  node;
   char*          name;
   char*          delim = " \t";    /* whitespace */

   if (fp == NULL) {
      fprintf(stderr, "ERROR: Unable to open file.\n" );
      exit(EXIT_FAILURE);
   }
   
   /* print index */
   fprintf(fp, "=== F_INDEX_NODE ===\n");
   fprintf(fp, "INDEX: %d\n", id);
   node = &(index->nodes[id]);
   name = strtok(node->name, delim);
   fprintf(fp, "%d\t%ld\t%s\t", id, node->offset, name );
   fprintf(fp, "\n");
}

/* *******************************************************************
 *    FUNC:    F_INDEX_Utest
 *    DESC:    F_INDEX unit test.
 * *******************************************************************/
void F_INDEX_Utest()
{
   F_INDEX* f_index = NULL;
   F_INDEX_Fasta_Build( f_index, "src/test.txt" );
   F_INDEX_Sort_by_Name(f_index);
   F_INDEX_Dump(f_index, stdout);

   for (int i = 0; i < f_index->N; i++) 
   {
      char* search_term = f_index->nodes[i].name;
      int   find        = F_INDEX_Search_Name(f_index, search_term);
      printf("INDEX: %d, FIND: %s\n", find, search_term);
   }
   exit(EXIT_SUCCESS);
}