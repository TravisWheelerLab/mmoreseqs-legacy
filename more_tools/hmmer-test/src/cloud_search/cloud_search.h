/*******************************************************************************
 *  FILE:      cloud_fwdback.c
 *  PURPOSE:   Runs Cloud Forward/Backward Algorithm.
 *             Cloud Search uses 3 stages:
 *             1.  Runs the "Cloud Search" Pruning in Forward/Backward directions,
 *                  which returns a reduced search space boundaries in the 
 *                  tradional dynamic programming matrix.
 *             2.  Forward/Backward passes are merged together.
 *             3.  Runs a bounded version of the traditional Forward/Backward.
 *
 *  BUG:       
 *******************************************************************************/

#ifndef _CLOUD_SEARCH_H_
#define _CLOUD_SEARCH_H_

/* === STAGE 1 === */
/* ------------------------------------------------------------------------------------------ *
 *  
 *  FUNCTION: cloud_Forward_Linear()
 *  SYNOPSIS: Perform Forward part of Cloud Search Algorithm.
 *            Traverses the dynamic programming matrix antidiagonally, running the
 *            Forward algorithm, starting at the Viterbi alignment beginning.  
 *            At the end of each antidiagonal, compares each cell against the current maximum
 *            scoring cell.  If cell falls below (MAX_SCORE * alpha), then cell is removed 
 *            from search space.  Terminates when reaches the final cell in dp matrix or 
 *            all cells in current antidiag have been pruned.  
 *            Stores final edgebound data in <edg>.
 *
 *  ARGS:      <query>     query sequence, 
 *             <target>    HMM model,
 *             <Q>         query length, 
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix (quadratic space),
 *             <st_MX3>    Normal State (Match, Insert, Delete) Matrix (linear space),
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix,
 *             <tr>        Traceback Data,
 *             <edg>       (OUTPUT) Edgebounds Tracker Data,
 *             <alpha>     Pruning Drop,
 *             <beta>      Number of Passes before Pruning Begins
 *
 *  RETURN:    No Return.
 *
 * ------------------------------------------------------------------------------------------- */
void cloud_Forward_Linear( const SEQUENCE*    query, 
                           const HMM_PROFILE* target,
                           const int          Q, 
                           const int          T, 
                           float*             st_MX, 
                           float*             st_MX3,
                           float*             sp_MX, 
                           const ALIGNMENT*   tr,
                           EDGEBOUNDS*        edg,
                           const float        alpha, 
                           const int          beta,
                           const bool         test);

/* ------------------------------------------------------------------------------------------ *
 *  
 *  FUNCTION: cloud_Backward_Linear()
 *  SYNOPSIS: Perform Backward part of Cloud Search Algorithm.
 *            Traverses the dynamic programming matrix antidiagonally, running the
 *            Forward algorithm, starting at the Viterbi alignment ending.  
 *            At the end of each antidiagonal, compares each cell against the current maximum
 *            scoring cell.  If cell falls below (MAX_SCORE * alpha), then cell is removed 
 *            from search space.  Terminates when reaches the final cell in dp matrix or 
 *            all cells in current antidiag have been pruned.  
 *            Stores final edgebound data in <edg>.
 *
 *  ARGS:      <query>     query sequence, 
 *             <target>    HMM model,
 *             <Q>         query length, 
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix (quadratic space),
 *             <st_MX3>    Normal State (Match, Insert, Delete) Matrix (linear space),
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix,
 *             <tr>        Traceback Data,
 *             <edg>       (OUTPUT) Edgebounds Tracker Data,
 *             <alpha>     Pruning Drop,
 *             <beta>      Number of Passes before Pruning Begins
 *
 *  RETURN:    No Return.
 *
 * ------------------------------------------------------------------------------------------- */
void cloud_Backward_Linear(const SEQUENCE*    query, 
                           const HMM_PROFILE* target,
                           const int          Q, 
                           const int          T, 
                           float*             st_MX, 
                           float*             st_MX3,
                           float*             sp_MX, 
                           const ALIGNMENT*   tr,
                           EDGEBOUNDS*        edg,
                           const float        alpha, 
                           const int          beta,
                           const bool         test);


/* === STAGE 2 === */
void EDGEBOUNDS_Reflect(EDGEBOUNDS *edg);

EDGEBOUNDS* EDGEBOUNDS_Merge(const int         Q, 
                             const int         T,
                             const EDGEBOUNDS* edg_1,
                             const EDGEBOUNDS* edg_2);

EDGEBOUNDS* EDGEBOUNDS_Reorient(const int         Q, 
                                const int         T,
                                const EDGEBOUNDS* edg_in);




/* === STAGE 3 === */
/* ------------------------------------------------------------------------------------------ *
 *
 *  FUNCTION: bound_Forward_Linear()
 *  SYNOPSIS: Perform Edge-Bounded Forward step of Cloud Search Algorithm.
 *            Runs traditional Forward-Backward Algorithm, but only performs
 *             computation on cells that fall within the bounds determined by
 *             the <edg> EDGEBOUNDS object, which stores a series of 
 *             (left-bound, right-bound) pairs sorted by row.  
 *            Normal state matrix is stored in linear space.
 *             <st_MX3> is size [3 * (Q + T + 1)]. Only requires size [2 * (T + 1)],
 *             but is reused from cloud_forward_().
 *            Final score produced by Forward is stored in <sc_final>.
 *
 *  ARGS:      <query>        Query sequence, 
 *             <target>       HMM model,
 *             <Q>            Query length, 
 *             <T>            Target length,
 *             <st_MX3>       Normal State (Match, Insert, Delete) Matrix (Linear Space),
 *             <st_MX>        Normal State (Match, Insert, Delete) Matrix (Quadratic Space),
 *             <sp_MX>        Special State (J,N,B,C,E) Matrix,
 *             <edg>          Edgebounds (stored row-wise)
 *             <test>         Toggles DEBUG output
 *             <sc_final>     (OUTPUT) Final Score 
 *
 *  RETURN:    Returns the final score of the Forward Algorithm.
 *
 * ------------------------------------------------------------------------------------------- */
float bound_Forward_Linear(const ESL_DSQ*     query, 
                           const P7_PROFILE*  target,
                           const int          Q, 
                           const int          T, 
                           float*             st_MX3,
                           float*             st_MX,
                           float*             sp_MX, 
                           const EDGEBOUNDS*  edg,
                           const bool         test,
                           float              *sc_final );

/* ------------------------------------------------------------------------------------------ *
 *
 *  FUNCTION: bound_Backward_Linear()
 *  SYNOPSIS: Perform Edge-Bounded Backward step of Cloud Search Algorithm.
 *            Runs traditional Backward Algorithm, but only performs
 *             computation on cells that fall within the bounds determined by
 *             the <edg> EDGEBOUNDS object, which stores a series of 
 *             (left-bound, right-bound) pairs sorted by row.  
 *            Normal state matrix is stored in linear space.
 *             <st_MX3> is size [3 * (Q + T + 1)]. Only requires size [2 * (T + 1)],
 *             but is reused from cloud_forward_Run3().
 *            Final score produced by Forward is stored in <sc_final>.
 *
 *  ARGS:      <query>        Query sequence, 
 *             <target>       HMM model,
 *             <Q>            Query length, 
 *             <T>            Target length,
 *             <st_MX3>       Normal State (Match, Insert, Delete) Matrix (Linear Space),
 *             <st_MX>        Normal State (Match, Insert, Delete) Matrix (Quadratic Space),
 *             <sp_MX>        Special State (J,N,B,C,E) Matrix,
 *             <edg>          Edgebounds (stored row-wise)
 *             <test>         Toggles DEBUG output
 *             <sc_final>     Final Score (OUTPUT)
 *
 *  RETURN:    Returns the final score of the Backward Algorithm.
 *
 * ------------------------------------------------------------------------------------------- */
float bound_Backward_Linear(  const ESL_DSQ*     query, 
                              const P7_PROFILE*  target,
                              const int          Q, 
                              const int          T, 
                              float*             st_MX3,
                              float*             st_MX,
                              float*             sp_MX, 
                              const EDGEBOUNDS*  edg,
                              const bool         test,
                              float              *sc_final );

#endif /* _CLOUD_SEARCH_H_ */
