/*******************************************************************************
 *  FILE:      	unordered_map.c
 *  PURPOSE:   	UNORDERED_MAP object.
 *  			when 
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
// #include "umap.h"

/* constructor */
UMAP* UMAP_Create( )
{
	UMAP* umap;
	int init_size = 16;

	umap = (UMAP*) malloc( sizeof(UMAP) );
	umap->N 		= 0;
	umap->Nalloc 	= init_size;
	umap->nodes		= (UMAP_NODE*) malloc( sizeof(UMAP_NODE) * init_size );

	return umap;
}

void UMAP_Add_Node( UMAP* 	umap, 
					char* 	key,
					DATA 	value )
{

}	

DATA UMAP_Get_Value( UMAP* 	umap,
					 char*	key )
{
	DATA data;
	data.i = 0;
	return data;
}
