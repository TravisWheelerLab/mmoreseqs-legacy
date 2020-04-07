/*******************************************************************************
 *  FILE:      structs.h
 *  PURPOSE:   Macros used by Cloud Search.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _STRUCTS_MACROS_H
#define _STRUCTS_MACROS_H

/* === MACRO FUNCTIONS === */

/* QUADRATIC DYNAMIC PROGRAMMING MATRIX - ACCESS MACROS */
/* generic access for any matrix */
#define ST_MX(mx,st,i,j)   ( mx[ ( st*(Q+1)*(T+1) ) + ( (i)*(T+1) ) + (j) ] )
/* match, insert, delete for st_MX matrix */
#define MMX(i,j)           ST_MX( st_MX, MAT_ST, i, j )
#define IMX(i,j)           ST_MX( st_MX, INS_ST, i, j )
#define DMX(i,j)           ST_MX( st_MX, DEL_ST, i, j )

/* LINEAR DYNAMIC PROGRAMMING MATRIX - ACCESS MACROS ( dim: 3 x (N+M) ) */
/* generic matrix */
#define ST_MX3(mx,st,i,j)  ( mx[ ( st*3*((T+1)+(Q+1)) ) + ( (i)*((T+1)+(Q+1)) ) + (j) ] )
/* match, insert, delete for st_MX3 matrix */
#define MMX3(i,j)          ST_MX3( st_MX3, MAT_ST, i, j )
#define IMX3(i,j)          ST_MX3( st_MX3, INS_ST, i, j )
#define DMX3(i,j)          ST_MX3( st_MX3, DEL_ST, i, j )

/* === MATRIX OBJECT FUNCTIONS === */
/* generic access for any matrix */
#define ST_MX_M(mx,st,i,j)  ( *MATRIX_3D_Get(mx,st,i,j) )
/* match, insert, delete for st_MX matrix */
#define MMX_M(i,j)           ST_MX_M( st_MX, MAT_ST, i, j )
#define IMX_M(i,j)           ST_MX_M( st_MX, INS_ST, i, j )
#define DMX_M(i,j)           ST_MX_M( st_MX, DEL_ST, i, j )

/* generic matrix */
#define ST_MX3_M(mx,st,i,j)  ( *MATRIX_3D_Get(mx,st,i,j) )
/* match, insert, delete for st_MX3 matrix */
#define MMX3_M(i,j)          ST_MX3_M( st_MX3, MAT_ST, i, j )
#define IMX3_M(i,j)          ST_MX3_M( st_MX3, INS_ST, i, j )
#define DMX3_M(i,j)          ST_MX3_M( st_MX3, DEL_ST, i, j )

#define SP_MX_M(mx,sp,i)     ( *MATRIX_2D_Get(mx,sp,i) )
#define XMX_M(sp,i)          SP_MX_M(sp_MX,sp,i)
/* ===================== */

/* SPECIAL STATE MATRIX MACROS */
#define SP_MX(mx,sp,i)     ( mx[ ((sp)*(Q+1)) + (i) ] )
#define XMX(sp,i)          SP_MX(sp_MX,sp,i)

/* TEST MATRIX */
#define TMX(i,j)           ( test_MX[ ((i)*(T+1)) + (j) ] ) 

/* TRANSITION SCORE, SPECIAL TRANSITION SCORE, MATCH SCORE, INSERT SCORE MACROS */
/* generic hmm profile functions */
#define TSC_HMM(prof,j,tr)    ( prof->hmm_model[j].trans[tr] )
#define XSC_HMM(prof,sp,tr)   ( prof->bg_model->spec[sp][tr] )
#define MSC_HMM(prof,j,A)     ( prof->hmm_model[j].match[A] )
#define ISC_HMM(prof,j,A)     ( prof->hmm_model[j].insert[A] )
/* target hmm profile functions */
#define TSC(j,tr)             ( target->hmm_model[j].trans[tr] )
#define XSC(sp,tr)            ( target->bg_model->spec[sp][tr] )
#define MSC(j,A)              ( target->hmm_model[j].match[A] )
#define ISC(j,A)              ( target->hmm_model[j].insert[A] )

/* DEBUG MACRO FOR RETREIVING VARIABLE NAME */
#define getName(var) #var

/* BASIC FUNCTION MACROS */
#define MAX(x,y)     (((x) > (y)) ? (x) : (y))
#define MIN(x,y)     (((x) < (y)) ? (x) : (y))
#define ABS(i)       (( (i) > (0) ? (i) : (-i) ))
/* check if two value are equal within tolerance */
#define CMP_TOL(i,j) (( fabs( (i) - (j) ) < tol ? 1 : 0 )) 

#define Test_IsLocal(mode)  (mode == MODE_MULTILOCAL || mode == MODE_UNILOCAL)
#define Test_IsMulti(mode)  (mode == MODE_MULTILOCAL || mode == MODE_MULTIGLOCAL)

#endif /* _STRUCTS_MACROS_H */
