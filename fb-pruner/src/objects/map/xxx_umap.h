/*******************************************************************************
 *  FILE:      	template_umap.h
 *  PURPOSE:   	XXX_UMAP object.
 *  					Unordered map of {key,value} pairs stored in hash table.
 * 					{key} is a STRING, {value} is type XXX.
 *
 *  AUTHOR:    	Dave Rich
 *******************************************************************************/

#ifndef XXX_UMAP_H
#define XXX_UMAP_H

/*    FUNCTION:   XXX_UMAP_Create()
 *    SYNOPSIS:   Create new XXX_UMAP struct.
 */
XXX_UMAP* XXX_UMAP_Create();

/*    FUNCTION:   XXX_UMAP_Destroy()
 *    SYNOPSIS:   Frees XXX_UMAP struct.
 */
XXX_UMAP* XXX_UMAP_Destroy();

/*    FUNCTION:   XXX_UMAP_Resize()
 *    SYNOPSIS:   Resize the XXX_UMAP struct.
 */
int XXX_UMAP_Resize( 	XXX_UMAP* 		umap, 
								int 				size );

/*    FUNCTION:   XXX_UMAP_Add()
 *    SYNOPSIS:   Add new {key,value} to {umap}.
 */
int XXX_UMAP_Add( 	XXX_UMAP* 		umap, 
							char* 			key,
							XXX*				value );

/*    FUNCTION:   XXX_UMAP_Get()
 *    SYNOPSIS:   Return {value} associated with {key,value} pair in {umap}.
 */
int XXX_UMAP_Get( 	XXX_UMAP* 	umap,
							char*			key,
							XXX* 			value );

#endif /* XXX_UMAP_H */