/*******************************************************************************
 *     FILE:   vizualizer.c
 *  PURPOSE:   Vizualizer 
 *
 *   AUTHOR:   Dave Rich
 *******************************************************************************/

#ifndef _VIZUALIZER_H
#define _VIZUALIZER_H

/*  FUNCTION:  VIZUALIZE_Edgebounds()
 *  SYNOPSIS:  
 */
void VIZUALIZER_Edgebounds(   EDGEBOUNDS*    edg,
                              char*          filename );

/*  FUNCTION:  VIZUALIZE_Matrix_2D()
 *  SYNOPSIS:  
 */
void VIZUALIZE_Matrix_2D(  MATRIX_2D*  mx,
                           char*       filename );

/*  FUNCTION:  VIZUALIZE_Matrix_3D()
 *  SYNOPSIS:  
 */
void VIZUALIZE_Matrix_3D(  MATRIX_3D*  mx,
                           char*       filename );

/*  FUNCTION:  VIZUALIZE_Matrix_3D()
 *  SYNOPSIS:  
 */
void VIZUALIZE_Matrix_3D_Sparse(    MATRIX_3D_SPARSE*    mx,
                                    char*                filename );

#endif /* _VIZUALIZER_H */