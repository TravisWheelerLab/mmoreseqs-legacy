/*******************************************************************************
 *  FILE:      structs.h
 *  PURPOSE:   Macros used by Cloud Search.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:       Lots.
 *******************************************************************************/

/* vectorization */
#include <xmmintrin.h>     /* SSE  */
#include <emmintrin.h>     /* SSE2 */

#ifndef _STRUCTS_MACROS_H
#define _STRUCTS_MACROS_H

#define TRUE  1
#define FALSE 0

/* === SET BUILD TYPE MACROS & FUNCTION COMPTILE-TIME OPTIONS  === */

/* set default debug build */
#ifndef DEBUG
#define DEBUG 	FALSE
#endif
/* set visualization build (subset of debug build) */
#ifndef VIZ
#define VIZ 	DEBUG
#endif
/* set memory checks (subset of debug build) */
#ifndef MEMCHECK
#define MEMCHECK 	DEBUG
#endif

/* whether to use function calls for matrix accesses or explicit array accesses */
#define MATRIX_FUNCTIONS 	FALSE

/* types of cloud pruning methods */
#define PRUNER_NONE  						0 
#define PRUNER_XDROP_EDGETRIM 				1
#define PRUNER_XDROP_BIFURCATE				2
#define PRUNER_DBL_XDROP_EDGETRIM_OR_DIE	3
/* if bifurcation is allowed, set max limit on forking paths */
#define MAX_BOUNDS_PER_ROW 		10
/* set default  of pruner method */
/* PRUNER METHODS: PRUNER_DBL_XDROP_EDGETRIM_OR_DIE, PRUNER_XDROP_EDGETRIM  */
#ifndef PRUNER
#define PRUNER  	PRUNER_DBL_XDROP_EDGETRIM_OR_DIE
#endif

/* types of cloud search: whether to store bounds as rows or antidiags */
#define CLOUD_NONE 		0
#define CLOUD_ROWS 		1
#define CLOUD_DIAGS 	2
/* set default cloud search method */
#ifndef CLOUD_METHOD
#define CLOUD_METHOD 	CLOUD_DIAGS
#endif

/* types of simd vectorization method */
#define SIMD_NONE 		0
#define SIMD_SSE		1
/* set default vectorization method */
#ifndef SIMD_METHOD
#define SIMD_METHOD 	SIMD_SSE
#endif

/* ============================================================================== */

/* values used for testing matrix accesses */
#define DIRTY_VAL 1.0
#define SCRUB_VAL 1.0

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
#define malloc_check(...) 	ERRORCHECK_malloc(__VA_ARGS__, ERRORCHECK_locate() )
#define realloc_check(...)	ERRORCHECK_realloc(__VA_ARGS__, ERRORCHECK_locate() )
#define fopen_check(...) 	ERRORCHECK_fopen(__VA_ARGS__, ERRORCHECK_locate() )
/* gets the location where error occurred */
#define ERRORCHECK_locate()			__FILE__, __LINE__, __FUNCTION__

/* === MATRIX FUNCTIONS AND MACROS === */

/* whether to access MATRIX_3D via function calls or direct data accesses */
#if ( MATRIX_FUNCTIONS == TRUE )
	/* generic access for MATRIX_3D via function call */
	#define MX_3D(mx, st, q_0, t_0)  	( *MATRIX_3D_Get( (mx), (st), (q_0), (t_0) ) )
#endif
#if ( MATRIX_FUNCTIONS == FALSE )
	/* generic access for MATRIX_3D via direct data access */
	#define MX_3D(mx, st, q_0, t_0)  	( mx->data[ ((st) * (mx->C * mx->N)) + ( (q_0) * (mx->N)) + (t_0) ] )
#endif

/* match, insert, delete for st_MX matrix (quadratic space matrix) */
#define MMX(q_0, t_0)           	MX_3D( st_MX, MAT_ST, (q_0), (t_0) )
#define IMX(q_0, t_0)           	MX_3D( st_MX, INS_ST, (q_0), (t_0) )
#define DMX(q_0, t_0)           	MX_3D( st_MX, DEL_ST, (q_0), (t_0) )
/* match, insert, delete for st_MX3 matrix (linear space matrix) */
#define MMX3(qx0, tx0)          	MX_3D( st_MX3, MAT_ST, (qx0), (tx0) )
#define IMX3(qx0, tx0)          	MX_3D( st_MX3, INS_ST, (qx0), (tx0) )
#define DMX3(qx0, tx0)          	MX_3D( st_MX3, DEL_ST, (qx0), (tx0) )

/* whether to use MATRIX_3D_SPARSE function calls or direct data accesses */
#if ( MATRIX_FUNCTIONS == TRUE )
	/* generic access for MATRIX_3D_SPARSE via function call */
	#define SMX(mx, st, qx0, tx0) 		( mx->data->data[ qx0 + (tx0 * NUM_NORMAL_STATES) + (st) ] )
#endif
#if ( MATRIX_FUNCTIONS == FALSE )
	/* generic access for MATRIX_3D_SPARSE via direct data access */
	#define SMX(mx, st, qx0, tx0) 		( mx->data->data[ (qx0) + ( (tx0) * NUM_NORMAL_STATES) + (st) ] )
#endif

/* match, insert, delete for st_SMX matrix (sparse matrix) */
/* NOTE: qx0 = bound mapped position, tx0 = t_0 - leftbound starting position */
#define MSMX(qx0, tx0) 			SMX( st_SMX, MAT_ST, (qx0), (tx0) )
#define ISMX(qx0, tx0) 			SMX( st_SMX, INS_ST, (qx0), (tx0) )
#define DSMX(qx0, tx0) 			SMX( st_SMX, DEL_ST, (qx0), (tx0) )

/* whether to access MATRIX_2D via function calls or direct data accesses */
#if ( MATRIX_FUNCTIONS == TRUE )
	/* generic access for MATRIX_2D via function call */
	#define MX_2D(mx, st, q_0)  	( *MATRIX_2D_Get( mx, st, q_0 ) )
#endif
#if ( MATRIX_FUNCTIONS == FALSE )
	/* generic access for MATRIX_2D via direct data access */
	#define MX_2D(mx, st, q_0)  	( mx->data[ ((st) * (mx->C)) + (q_0) ] )
#endif

/* special state */
#define XMX(st, q_0)          	MX_2D( sp_MX, st, q_0)

/* === TRANSITION SCORE, SPECIAL TRANSITION SCORE, MATCH SCORE, INSERT SCORE MACROS === */
/* generic hmm profile functions */
#define TSC_HMM(prof, j, tr)    ( prof->hmm_model[j].trans[tr] )
#define XSC_HMM(prof, sp, tr)   ( prof->bg_model->spec[sp][tr] )
#define MSC_HMM(prof, j, A)     ( prof->hmm_model[j].match[A] )
#define ISC_HMM(prof, j, A)     ( prof->hmm_model[j].insert[A] )
/* target hmm profile functions */
#define TSC(j, tr)             	( target->hmm_model[j].trans[tr] )
#define XSC(sp, tr)            	( target->bg_model->spec[sp][tr] )
#define MSC(j, A)              	( target->hmm_model[j].match[A] )
#define ISC(j, A)              	( target->hmm_model[j].insert[A] )

/* edgebounds access */
#define EDG_X(edg, i) 	 	( *EDGEBOUNDS_Get( edg, i ) )

/* logarithmic sum */
#define LOGSUM(a, b) 		( logsum( a, b ) )

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

/* === SIMD VECTORIZATION FUNCTIONS === */
#if ( SIMD_METHOD == SIMD_SSE )
	#define VECTOR_WIDTH 	128
	typedef __m128i 		__VECTOR_INT;
	typedef __m128 			__VECTOR_FLT;	
	typedef __m128d			__VECTOR_DBL;
#endif

#define BITS_PER_BYTE		8
#define CHAR_BITS			1 * BITS_PER_BYTE
#define INT_BITS 			4 * BITS_PER_BYTE
#define FLOAT_BITS 			4 * BITS_PER_BYTE

#define CHARS_PER_VEC 	 	VECTOR_WIDTH / CHAR_BITS
#define INTS_PER_VEC 		VECTOR_WIDTH / INT_BITS
#define FLOATS_PER_VEC 		VECTOR_WIDTH / FLOAT_BITS

/* number of float vectors per sequence (striped implementation requires at least two vectors) */
#define NUM_VEC_FOR_SEQ(seq_length, data_per_vec)  calc_Max( 2, (seq_length / data_per_vec) + 1 )

/* matrix access functions */
#define MX_3D_VEC(mx,st,i,j)	
#define MMV(i,j) 				MX_3D_VEC(  )

#endif /* _STRUCTS_MACROS_H */
