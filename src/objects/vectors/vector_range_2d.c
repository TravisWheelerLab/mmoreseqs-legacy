/*******************************************************************************
 *  @file vector_range_2d.c
 *  @brief 2D RANGE VECTOR objects
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

/* local imports */
#include "structs.h"
#include "utilities.h"
#include "objects.h"

/* header */
#include "vector_range_2d.h"


/* constructor */
VECTOR_RANGE_2D* VECTOR_RANGE_2D_Create()
{
   VECTOR_RANGE_2D *vec;
   const int init_size = 8;
   vec = (VECTOR_RANGE_2D *) malloc( sizeof(VECTOR_RANGE_2D) );
   vec->id = VECTOR_INT_Create();
   vec->data = (VECTOR_RANGE *) malloc( sizeof(VECTOR_RANGE) * init_size );
   vec->N = 0;
   vec->Nalloc = init_size;
   return vec;
}

/* constructor with initial size and fill  */
VECTOR_RANGE_2D* VECTOR_RANGE_2D_Create_Size( const int init_size )
{
   VECTOR_RANGE_2D *vec;
   vec = (VECTOR_RANGE_2D *) malloc( sizeof(VECTOR_RANGE_2D) );
   vec->id = VECTOR_INT_Create();
   vec->data = (VECTOR_RANGE *) malloc( sizeof(VECTOR_RANGE) * (init_size + 1) );
   for (int i = 0; i < init_size; i++) {
      VECTOR_RANGE *subvec = VECTOR_RANGE_Create();
      VECTOR_RANGE_2D_Pushback( vec, *subvec );
   }
   vec->N = init_size;
   vec->Nalloc = init_size + 1;
   return vec;
}

/* destructor */
void VECTOR_RANGE_2D_Destroy( VECTOR_RANGE_2D *vec )
{
   VECTOR_INT_Destroy(vec->id);
   VECTOR_RANGE_Destroy(vec->data);
   free(vec);
}

/* deep copy */
VECTOR_RANGE_2D* VECTOR_RANGE_2D_Copy( VECTOR_RANGE_2D *src )
{
   VECTOR_RANGE_2D *vec;
   vec = (VECTOR_RANGE_2D *) malloc( sizeof(VECTOR_RANGE_2D) );
   /* copy base data */
   memcpy( vec, src, sizeof(VECTOR_RANGE_2D) );
   /* copy variable-sized data */
   vec->data = (VECTOR_RANGE *) malloc( sizeof(VECTOR_RANGE) * src->Nalloc );
   memcpy( vec->data, src->data, sizeof(VECTOR_RANGE) * src->N );

   for (int i = 0; i < vec->N; i++) {
      vec->data[i].data = (RANGE *) malloc( sizeof(RANGE) * src->data[i].Nalloc );
      memcpy( vec->data[i].data, src->data[i].data, sizeof(RANGE) * src->data[i].N );
   }

   return vec;
}

/* resize the array */
void VECTOR_RANGE_2D_Resize( VECTOR_RANGE_2D *vec, float growth_factor )
{
   vec->data = (VECTOR_RANGE *) realloc( vec->data, sizeof(VECTOR_RANGE) * vec->Nalloc * growth_factor );
   vec->Nalloc *= growth_factor;
}

/* push element onto end of array */
void VECTOR_RANGE_2D_Pushback( VECTOR_RANGE_2D *vec, VECTOR_RANGE val )
{
   /* if array is full, resize */
   if (vec->N >= vec->Nalloc - 1) {
      VECTOR_RANGE_2D_Resize( vec, 2 );
   }

   vec->data[vec->N] = val;
   vec->N++;
}

/* pop element from end of array */
VECTOR_RANGE VECTOR_RANGE_2D_Pop( VECTOR_RANGE_2D *vec )
{
   VECTOR_RANGE tmp = vec->data[vec->N-1];
   vec->N -= 1;

   /* if array is less than half used, resize */
   if (vec->N < vec->Nalloc / 2) {
      VECTOR_RANGE_2D_Resize( vec, 0.5 );
   }

   return tmp;
}

/* set data at index (no bound checks) */
void VECTOR_RANGE_2D_Set( VECTOR_RANGE_2D *vec, int idx, VECTOR_RANGE val )
{
   vec->data[idx] = val;
}

/* get data at index (no checks) */
VECTOR_RANGE VECTOR_RANGE_2D_Get( VECTOR_RANGE_2D *vec, int idx )
{
   return vec->data[idx];
}

/* compare two VECTOR_RANGE_2D objects */
int VECTOR_RANGE_2D_Compare( VECTOR_RANGE_2D *vecA, VECTOR_RANGE_2D *vecB )
{
   // for (int i = 0; i < vecA->N; i++) {
   //    if ( vecA->data[i] != vecB->data[i] ) {
   //       if ( vecA->data[i] > vecB->data[i] ) {
   //          return 1;
   //       } else {
   //          return -1;
   //       }
   //    }
   // }
   return 0;
}

/* merge current diagonal bound into vectors in the forward direction */
void VECTOR_RANGE_2D_MergeFwd( VECTOR_RANGE_2D *vec, const BOUND bnd ) 
{
   VECTOR_RANGE row;
   const int tol = 2;
   int d = bnd.id;
   /* find where diagonal intersects each row */
   for (int i = bnd.lb; i < bnd.rb; i++) {
      int j = d - i;
      row = vec->data[i];
      /* check if cell is adjacent to current open range */
      if ( j - row.data[(row.N)-1].end < tol ) 
      {
         row.data[(row.N)-1].end = j + 1;
      }
      /* otherwise, create a new open range on row */
      else
      {
         VECTOR_RANGE_Pushback( &row, (RANGE){j,j+1} );
      }
   }
}

/* convert 2d vector of flat list of edgebounds with head pointers */
EDGEBOUNDS* VECTOR_RANGE_2D_Convert_to_Edgebound( VECTOR_RANGE_2D *vec ) 
{
//    EDGEBOUNDS   *edg = EDGEBOUNDS_Create();
//    VECTOR_RANGE row;
//    RANGE        range;

//    for (int i = 0; i < vec->N; i++) 
//    {
//       row = vec->data[i];
//       EDGEBOUNDS_Pushback_Head( edg, i, edg->N );
//       for (int j = 0; j < row.N; j++)
//       {
//          range = row.data[j];
//          EDGEBOUNDS_Pushback( edg, (BOUND){i, range.beg, range.end} );
//       }
//    }

//    return edg;
}

void VECTOR_RANGE_2D_Dump(FILE *fp, VECTOR_RANGE_2D *vec)
{
   printf("=== VECTOR RANGE 2D ===\n");
   for (int i = 0; i < vec->N; i++) 
   {
      VECTOR_RANGE subvec = vec->data[i];
      printf("[%d]=>%d:\n", i, subvec.N);

      for (int j = 0; j < subvec.N; j++)
      {
         RANGE r = subvec.data[j];
         printf("\t[%d] (%d,%d)\n", j, r.beg, r.end);
      }
      printf("\n");
   }
   printf("\n");
}


void VECTOR_RANGE_2D_UnitTest()
{

   printf("Begin unit test...\n");
   srand(time(0));
   VECTOR_RANGE_2D *vec = VECTOR_RANGE_2D_Create();
   printf("Size: %d/%d\n", vec->N, vec->Nalloc);

   printf("Loading 2-D Vector...\n");
   for (int i = 0; i < 20; i++) 
   {
      VECTOR_RANGE *subvec = VECTOR_RANGE_Create();
      int num_range = ( rand() % 10 ) + 10;
      for (int j = 0; j < num_range; j++) {
         int a = ( rand() % (1000-100) ) + 100;
         int b = ( rand() % (1000-a) ) + a;
         RANGE r = (RANGE){a,b};
         VECTOR_RANGE_Pushback(subvec, r);
      }
      VECTOR_RANGE_2D_Pushback(vec, *subvec);
   }
   VECTOR_RANGE_2D_Dump(stdout, vec);
   exit(EXIT_SUCCESS);
}

