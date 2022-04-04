/*******************************************************************************
 *  - FILE:  work_scoring.h
 *  - DESC:  Workflow Subroutines for Computing Final Scores.
 *******************************************************************************/

#ifndef _WORK_SCORING_H
#define _WORK_SCORING_H

/*! FUNCTION:  WORK_construct_scores()
 *  SYNOPSIS:  Construct final scores and evalues from results.
 */
void WORK_construct_scores(WORKER* worker);

/*! FUNCTION:  WORK_construct_scores_bydom()
 *  SYNOPSIS:  Construct final scores and evalues from domain results.
 */
void WORK_construct_scores_bydom(WORKER* worker);

#endif /* _WORK_SCORING_H */
