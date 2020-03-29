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

/* MATCH, INSERT, DELETE, SPECIAL DP MATRIX ACCESS MACROS */
#define ST_MX(mx,st,i,j)   ( mx[ ( st*(Q+1)*(T+1) ) + ( (i)*(T+1) ) + (j) ] )
#define MMX(i,j)           ST_MX( st_MX, MAT_ST, i, j)
#define IMX(i,j)           ST_MX( st_MX, INS_ST, i, j)
#define DMX(i,j)           ST_MX( st_MX, DEL_ST, i, j)
/* MATCH, INSERT, DELETE, SPECIAL ANTI-DIAG ACCESS MACROS (d = diag, i = offset) */
#define ST_MX3(mx,st,i,j)  ( mx[ ( st*3*((T+1)+(Q+1)) ) + ( (i)*((T+1)+(Q+1)) ) + (j) ] )
#define MMX3(i,j)          ST_MX3( st_MX3, MAT_ST, i, j)
#define IMX3(i,j)          ST_MX3( st_MX3, INS_ST, i, j)
#define DMX3(i,j)          ST_MX3( st_MX3, DEL_ST, i, j)

/* SPECIAL STATE MATRIX MACROS */
#define SP_MX(mx,sp,i)     (mx[ ((sp)*(Q+1)) + (i) ])
#define XMX(sp,i)          (sp_MX[ ((sp)*(Q+1)) + (i) ])

/* TEST MATRIX */
#define TMX(i,j)           (test_MX[ ((i)*(T+1)) + (j) ])

/* TRANSITION SCORE, SPECIAL TRANSITION SCORE, MATCH SCORE, INSERT SCORE MACROS */
#define TSC(j,tr)          (target->hmm_model[j].trans[tr])
// #define TSC(j,tr)          (target->bg_model->trans[tr])
#define XSC(sp,tr)         (target->bg_model->spec[sp][tr])
#define MSC(j,A)           (target->hmm_model[j].match[A])
#define ISC(j,A)           (target->hmm_model[j].insert[A])

/* DEBUG MACRO FOR RETREIVING VARIABLE NAME */
#define getName(var) #var

/* BASIC FUNCTION MACROS */
#define MAX(x,y)     (((x) > (y)) ? (x) : (y))
#define MIN(x,y)     (((x) < (y)) ? (x) : (y))
#define ABS(i)       (( (i) > (0) ? (i) : (-i) ))
/* check if two value are equal within tolerance */
#define CMP_TOL(i,j) (( fabs( (i) - (j) ) < tol ? 1 : 0 )) 

/* DEFINED CONSTANTS */
#define CONST_LOG2 0.69314718055994529
#define SCALE_FACTOR 1000
#define INF INFINITY
#define INT_MIN -2147483648

#define Test_IsLocal(mode)  (mode == MODE_MULTILOCAL || mode == MODE_UNILOCAL)
#define Test_IsMulti(mode)  (mode == MODE_MULTILOCAL || mode == MODE_MULTIGLOCAL)

#endif /* _STRUCTS_MACROS_H */
