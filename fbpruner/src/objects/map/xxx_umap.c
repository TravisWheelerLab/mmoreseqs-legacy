/*******************************************************************************
 *  FILE:      	template_umap.c
 *  PURPOSE:   	XXX_UMAP object.
 *  					Unordered map of {key,value} pairs stored in hash table.
 * 					{key} is a STRING, {value} is type XXX.
 *
 *  AUTHOR:    	Dave Rich
 *  BUG:
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

/* local imports */
#include "structs.h"
#include "utilities.h"
#include "objects.h"

/* header */
#include "xxx_umap.h"

# define HASH_TABLE_SIZE_LIMIT 

/*    FUNCTION:   XXX_UMAP_Create()
 *    SYNOPSIS:   Create new XXX_UMAP struct.
 */
XXX_UMAP* XXX_UMAP_Create()
{
	XXX_UMAP* umap;
	int init_size = 13;

	umap 				= (XXX_UMAP*) ERROR_malloc( sizeof(XXX_UMAP) );
	umap->N 			= 0;
	umap->Nalloc 	= init_size;
	umap->nodes		= (XXX_UMAP_NODE*) ERROR_malloc( sizeof(XXX_UMAP_NODE) * init_size );

	return umap;
}

/*    FUNCTION:   XXX_UMAP_Destroy()
 *    SYNOPSIS:   Frees XXX_UMAP struct.
 */
XXX_UMAP* XXX_UMAP_Destroy( XXX_UMAP* umap )
{
	ERROR_free( umap->nodes );
	ERROR_free( umap );
}

/*    FUNCTION:   XXX_UMAP_Resize()
 *    SYNOPSIS:   Resize the XXX_UMAP struct.
 */
int XXX_UMAP_Resize( 	XXX_UMAP* 		umap, 
								int 				size )
{

}

/*    FUNCTION:   XXX_UMAP_Resize()
 *    SYNOPSIS:   Resize the XXX_UMAP struct.
 */
int XXX_UMAP_GrowTo( 	XXX_UMAP* 		umap, 
								int 				size )
{
	if ( umap->N < size ) XXX_UMAP_Resize( umap, size );
}

/*    FUNCTION:   XXX_UMAP_Add()
 *    SYNOPSIS:   Add new {key,value} to hash table.
 */
int XXX_UMAP_Add( 	XXX_UMAP* 		umap, 
							char* 			key,
							XXX* 				value )
{

}	

/*    FUNCTION:   XXX_UMAP_Delete()
 *    SYNOPSIS:   Add {key} from Hash Table.  Returns 
 */
int XXX_UMAP_Remove( 	XXX_UMAP* 		umap, 
							 	char* 			key )
{
	
}

/*    FUNCTION:   XXX_UMAP_Add()
 *    SYNOPSIS:   Add new {key,value} to hash table.
 */
int XXX_UMAP_Get( 	XXX_UMAP* 	umap,
							char*			key,
							XXX* 			value )
{
	XXX data;
	return data;
}

/*    FUNCTION:   XXX_UMAP_Hash_Function()
 *    SYNOPSIS:   
 */
void XXX_UMAP_Hash_Function( 	XXX_UMAP* 		umap, 
										char* 			key )
{

}	
