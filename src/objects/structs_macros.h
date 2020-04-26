/*******************************************************************************
 *  FILE:      structs.h
 *  PURPOSE:   Macros used by Cloud Search.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

#ifndef _STRUCTS_MACROS_H
#define _STRUCTS_MACROS_H

/* === BUILD TYPE MACROS === */
/* set debug macros */
#ifndef DEBUG
#define DEBUG true
#endif

/* debug print (eliminated from code when not debugging) */
#if DEBUG
#define DBG_PRINTF(...) 	printf(__VA_ARGS__)
#else
#define DBG_PRINTF(...) 
#endif
#if DEBUG
#define DBG_FPRINTF(...) 	fprintf(__VA_ARGS__)
#else
#define DBG_FPRINTF(...) 
#endif

/* determines whether output is printed based on verbose_level */
#define printf_v(v_min,...) 	if ( args->verbose_level >= v_min ) printf(__VA_ARGS__)
#define printf_vnone(...) 		printf_v(VERBOSE_NONE, __VA_ARGS__)
#define printf_vlo(...) 		printf_v(VERBOSE_LOW,  __VA_ARGS__)
#define printf_vhi(...) 		printf_v(VERBOSE_HIGH, __VA_ARGS__)
#define printf_vall(...) 		printf_v(VERBOSE_ALL,  __VA_ARGS__)

/* errorchecking macros */
#define malloc_check(...) 	ERRORCHECK_malloc(__VA_ARGS__, __FILE__, __LINE__, __FUNCTION__)
#define realloc_check(...)	ERRORCHECK_realloc(__VA_ARGS__, __FILE__, __LINE__, __FUNCTION__)
#define fopen_check(...) 	ERRORCHECK_fopen(__VA_ARGS__, __FILE__, __LINE__, __FUNCTION__)
/* gets the location where error occurred */
#define LOCATION 			__FILE__, __LINE__, __FUNCTION__

/* === MATRIX FUNCTIONS === */
/* generic access for any 3d matrix */
#define MX_3D(mx,st,i,j)  	( *MATRIX_3D_Get( mx, st, i, j ) )
/* match, insert, delete for st_MX matrix */
#define MMX(i,j)           	MX_3D( st_MX, MAT_ST, i, j )
#define IMX(i,j)           	MX_3D( st_MX, INS_ST, i, j )
#define DMX(i,j)           	MX_3D( st_MX, DEL_ST, i, j )
/* match, insert, delete for st_MX3 matrix */
#define MMX3(i,j)          MX_3D( st_MX3, MAT_ST, i, j )
#define IMX3(i,j)          MX_3D( st_MX3, INS_ST, i, j )
#define DMX3(i,j)          MX_3D( st_MX3, DEL_ST, i, j )

/* generic access for any 2d matrix */
#define MX_2D(mx,sp,i)     ( *MATRIX_2D_Get(mx,sp,i) )
#define XMX(sp,i)          MX_2D(sp_MX,sp,i)

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

/* edgebounds access */
#define EDG_X(edg,i) 	 	( *EDGEBOUNDS_Get( edg, i ) )

/* logarithmic sum */
#define LOGSUM(a,b) 		( logsum( a, b ) )

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
