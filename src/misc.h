/******************************************************************************
 *  @file misc.h
 *  @brief Miscellaneous helper functions.
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/

#ifndef _MISC_H
#define _MISC_H

// constants
#define LOGSUM_SCALE 1000.f
#define LOGSUM_TBL 16000

/* Min/Max Fns for Viterbi */
// static inline
float calc_Max (float x, float y);
// static inline
float calc_Min (float x, float y);

/* Logsum fns for Forward/Backward */
void init_Logsum ();
// static inline
float calc_Logsum (float x, float y);
// static inline
float calc_Logsum_exact (float x, float y);
void print_Logsum ();
int get_str_len(char* str) ;

/* DP Matrix fns */
void dp_matrix_Print (const int Q, const int T, 
                      const float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                      const float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ]);
void dp_matrix_Clear (const int Q, const int T, 
                      float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                      float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ]);
void dp_matrix_Clear3 (const int Q, const int T,
                      float st_MX3[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                      float sp_MX[ NUM_SPECIAL_STATES * (Q + 1) ]);
void dp_matrix_Clear_X (const int Q, const int T, 
                       float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                       float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ],
                       float val);
void dp_matrix_Save (const int Q, const int T, 
                     const float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                     const float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ],
                     const char *_filename_);
void dp_matrix_trace_Save (const int Q, const int T, 
                        const float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                        const float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ],
                        const TRACEBACK *tr,
                        const char *_filename_);
void test_matrix_Print (const int Q, const int T, 
                        const float test_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ]);
int dp_matrix_Compare (const int Q, const int T,
                        float st_MX_1[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                        float sp_MX_1[ NUM_SPECIAL_STATES * (Q + 1) ],
                        float st_MX_2[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                        float sp_MX_2[ NUM_SPECIAL_STATES * (Q + 1) ] );
void dp_matrix_Copy (const int Q, const int T,
                     float st_MX_src[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                     float sp_MX_src[ NUM_SPECIAL_STATES * (Q + 1) ],
                     float st_MX_dest[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                     float sp_MX_dest[ NUM_SPECIAL_STATES * (Q + 1) ] );

#endif /* _MISC_H */