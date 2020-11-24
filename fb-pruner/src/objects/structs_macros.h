/*******************************************************************************
 *  FILE:      structs.h
 *  PURPOSE:   Macros used by fb-pruner.  Some set by Makefile.
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

#define MACRO_XSTR(val) 	MACRO_STR(val)
#define MACRO_STR(val) 		#val

/* used for storing potential string paths. Max determined by Linux */
#ifndef MAX_PATH_LEN
#define MAX_PATH_LEN 		4096
#endif 

/* === SET BUILD TYPE MACROS & FUNCTION COMPILE-TIME OPTIONS  === */

/* === VERSION === */
#define BUILD_PROGRAM      "MMORE-SEQS // MMSEQS-PLUS // FB-PRUNER"
#define BUILD_VERSION      "0.1"
#define BUILD_NAME         "tbd"
#define BUILD_DATE         "Aug 2020"
#define BUILD_DESCRIPT     "Heuristic Pruned Forward-Backward for Faster Profile-to-Sequence Search"

/* === INSTALL LOCATION === */
#ifndef PROJECT_LOCATION
#define PROJECT_LOCATION   ./
#endif

#ifndef INSTALL_LOCATION 
#define INSTALL_LOCATION	./
#endif

/* set default debug build */
#ifndef DEBUG
#define DEBUG 		FALSE
#endif
/* set visualization build (subset of debug build) */
#ifndef VIZ
#define VIZ 		DEBUG
#endif
/* set memory checks (subset of debug build) */
#ifndef MEMCHECK
#define MEMCHECK 	DEBUG
#endif

/* whether to use function calls for matrix accesses or explicit array accesses */
#define MATRIX_FUNCTIONS 	FALSE

/* types of cloud pruning methods */
#define PRUNER_NONE  							0 
#define PRUNER_XDROP_EDGETRIM 					1
#define PRUNER_XDROP_BIFURCATE					2
#define PRUNER_DBL_XDROP_EDGETRIM_OR_DIE		3
/* if bifurcation is allowed, set max limit on forking paths */
#ifndef MAX_BOUNDS_PER_ROW
#define MAX_BOUNDS_PER_ROW 						10
#endif
/* set default  of pruner method */
/* PRUNER METHODS: PRUNER_DBL_XDROP_EDGETRIM_OR_DIE, PRUNER_XDROP_EDGETRIM  */
#ifndef PRUNER
#define PRUNER  	PRUNER_DBL_XDROP_EDGETRIM_OR_DIE
#endif

/* types of cloud search: whether to store bounds as rows or antidiags */
#define CLOUD_NONE 			0
#define CLOUD_ROWS 			1
#define CLOUD_DIAGS 		2
/* set default cloud search method */
#ifndef CLOUD_METHOD
#define CLOUD_METHOD 		CLOUD_DIAGS
#endif

/* types of simd vectorization method */
#define SIMD_NONE 			0
#define SIMD_SSE			1
/* set default vectorization method */
#ifndef SIMD_METHOD
#define SIMD_METHOD 	SIMD_SSE
#endif

/* ============================================================================== */

/* values used for testing matrix accesses */
#define DIRTY_VAL 	1.0
#define SCRUB_VAL 	1.0

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
#define printf_v(v_threshold,...) 	if ( args->verbose_level >= v_threshold ) printf(__VA_ARGS__)
#define printf_vnone(...) 			printf_v(VERBOSE_NONE, __VA_ARGS__)
#define printf_vlo(...) 			printf_v(VERBOSE_LOW,  __VA_ARGS__)
#define printf_vhi(...) 			printf_v(VERBOSE_HIGH, __VA_ARGS__)
#define printf_vall(...) 			printf_v(VERBOSE_ALL,  __VA_ARGS__)

/* determines whether output is printed based on verbose_level */
#define fprintf_v(v_threshold,...) 	if ( args->verbose_level >= v_threshold ) fprintf(__VA_ARGS__)
#define fprintf_vnone(...) 		fprintf_v(VERBOSE_NONE, __VA_ARGS__)
#define fprintf_vlo(...) 			fprintf_v(VERBOSE_LOW,  __VA_ARGS__)
#define fprintf_vhi(...) 			fprintf_v(VERBOSE_HIGH, __VA_ARGS__)
#define fprintf_vall(...) 			fprintf_v(VERBOSE_ALL,  __VA_ARGS__)

/* error-checking macros */
#define ERROR_alloc(...) 			ERRORCHECK_alloc(__VA_ARGS__, ERRORCHECK_locate() )
#define ERROR_malloc(...) 			ERRORCHECK_malloc(__VA_ARGS__, ERRORCHECK_locate() )
#define ERROR_realloc(...) 			ERRORCHECK_realloc(__VA_ARGS__, ERRORCHECK_locate() )
#define ERROR_free(...) 			ERRORCHECK_free(__VA_ARGS__, ERRORCHECK_locate() )
#define ERROR_fopen(...) 			ERRORCHECK_fopen(__VA_ARGS__, ERRORCHECK_locate() )
#define ERROR_fclose(...) 			ERRORCHECK_fclose(__VA_ARGS__, ERRORCHECK_locate() )

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

/* match, insert, delete for st_MX matrix (quadratic space matrix) (specify matrix) */
#define MMX_X(mx, q_0, t_0)      MX_3D( mx, MAT_ST, (q_0), (t_0) )
#define IMX_X(mx, q_0, t_0)      MX_3D( mx, INS_ST, (q_0), (t_0) )
#define DMX_X(mx, q_0, t_0)      MX_3D( mx, DEL_ST, (q_0), (t_0) )
/* match, insert, delete for st_MX3 matrix (linear space matrix) (specify matrix) */
#define MMX3_X(mx, qx0, tx0)     MX_3D( mx, MAT_ST, (qx0), (tx0) )
#define IMX3_X(mx, qx0, tx0)     MX_3D( mx, INS_ST, (qx0), (tx0) )
#define DMX3_X(mx, qx0, tx0)     MX_3D( mx, DEL_ST, (qx0), (tx0) )

/* whether to use MATRIX_3D_SPARSE function calls or direct data accesses */
#if ( MATRIX_FUNCTIONS == TRUE )
	/* generic access for MATRIX_3D_SPARSE via function call */
	#define SMX(mx, st, qx0, tx0) 		( mx->data->data[ (qx0) + ( (tx0) * NUM_NORMAL_STATES ) + (st) ] )
#endif
#if ( MATRIX_FUNCTIONS == FALSE )
	/* generic access for MATRIX_3D_SPARSE via direct data access */
	#define SMX(mx, st, qx0, tx0) 		( mx->data->data[ (qx0) + ( (tx0) * NUM_NORMAL_STATES ) + (st) ] )
#endif

/* match, insert, delete for st_SMX matrix (sparse matrix) */
/* NOTE: qx0 = bound mapped position, tx0 = t_0 - leftbound starting position */
#define MSMX(qx0, tx0) 			SMX( st_SMX, MAT_ST, (qx0), (tx0) )
#define ISMX(qx0, tx0) 			SMX( st_SMX, INS_ST, (qx0), (tx0) )
#define DSMX(qx0, tx0) 			SMX( st_SMX, DEL_ST, (qx0), (tx0) )
/* match, insert, delete for st_SMX matrix (sparse matrix) */
#define MSMX_X(mx, qx0, tx0) 	SMX( mx, MAT_ST, (qx0), (tx0) )
#define ISMX_X(mx, qx0, tx0) 	SMX( mx, INS_ST, (qx0), (tx0) )
#define DSMX_X(mx, qx0, tx0) 	SMX( mx, DEL_ST, (qx0), (tx0) )

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
#define XMX_X(mx, st, q_0)       MX_2D( mx, st, q_0)

/* === TRANSITION SCORE, SPECIAL TRANSITION SCORE, MATCH SCORE, INSERT SCORE MACROS === */
/* target hmm profile functions */
#define TSC(t_0, tr)            	( target->hmm_model[t_0].trans[tr] )
#define XSC(sp, tr)            	( target->bg_model->spec[sp][tr] )
#define MSC(t_0, A)             	( target->hmm_model[t_0].match[A] )
#define ISC(t_0, A)            	( target->hmm_model[t_0].insert[A] )
/* generic hmm profile functions */
#define TSC_X(prof, t_0, tr)   	( prof->hmm_model[t_0].trans[tr] )
#define XSC_X(prof, sp, tr)    	( prof->bg_model->spec[sp][tr] )
#define MSC_X(prof, t_0, A)    	( prof->hmm_model[t_0].match[A] )
#define ISC_X(prof, t_0, A)    	( prof->hmm_model[t_0].insert[A] )

/* edgebounds access */
#define EDG_X(edg, i) 	 	( *EDGEBOUNDS_Get( edg, i ) )

/* vector access (vector, index) */
#define VEC_X(vec, i) 		( vec->data[i] )

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
	typedef __m128i 			__VECTOR_INT;
	typedef __m128 			__VECTOR_FLT;	
	typedef __m128d			__VECTOR_DBL;
#endif

#define BITS_PER_BYTE		8
#define CHAR_BITS				1 * BITS_PER_BYTE
#define INT_BITS 				4 * BITS_PER_BYTE
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
