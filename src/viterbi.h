/*******************************************************************************
 *  @file viterbi.h
 *  @brief The Viterbi Algorithm for Sequence Alignment Search.
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/

#ifndef _VITERBI_H
#define _VITERBI_H

float viterbi_Run (const SEQ* query, 
                  const HMM_PROFILE* target, 
                  int Q, int T, 
                  float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                  float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                  RESULTS* res,
                  TRACEBACK* tr);

void traceback_Build (const SEQ* query, 
                     const HMM_PROFILE* target, 
                     int Q, int T, 
                     float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                     float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ],
                     TRACEBACK* tr);

void traceback_Append (TRACEBACK* tr,
                     int st,
                     int i,
                     int j);

void traceback_Reverse (TRACEBACK* tr);

void traceback_Show (const int Q, const int T, 
                     float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                     float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ],
                     TRACEBACK *tr);

void traceback_Print (TRACEBACK *tr);

void traceback_Save(TRACEBACK *tr,
                    const char *_filename_);

#endif /* _VITERBI_H */