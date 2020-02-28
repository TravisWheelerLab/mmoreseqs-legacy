/*******************************************************************************
 *  @file cloud_search.h
 *  @brief Testing for navigating through the matrices.
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/


#ifndef _TESTING_H
#define _TESTING_H

/* === INCLUDES === */
// #include "objects/structs.h"
// #include "objects/edgebound.h"
// #include "objects/alignment.h"


/* === FUNCTIONS === */
void fwd_test_cycle(int Q, int T,
                   float *st_MX,
                   float *sp_MX,
                   ALIGNMENT *aln);

void bck_test_cycle(int Q, int T,
                    float *st_MX,
                    float *sp_MX,
                    ALIGNMENT *aln);

void fwd_test_cycle3(int Q, int T, 
                     float* st_MX, 
                     float* st_MX3,
                     float* sp_MX, 
                     ALIGNMENT *aln );

void bck_test_cycle3(int Q, int T, 
                     float* st_MX, 
                     float* st_MX3,
                     float* sp_MX, 
                     ALIGNMENT *aln );

int cloud_Fill(int Q, int T,
               float *st_MX,
               float *sp_MX,
               EDGEBOUNDS* edg,
               float val, 
               int mode);

int cloud_Solid_Fill(int Q, int T,
                     float *st_MX,
                     float *sp_MX,
                     EDGEBOUNDS* edg,
                     float val, 
                     int mode);

int cloud_Cell_Count(int Q, int T,
                     float *st_MX,
                     float *sp_MX);

/* DP MATRIX FUNCTIONS */
void dp_matrix_Print (const int Q, const int T, 
                      const float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                      const float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ]);
void dp_matrix_Print3 (const int Q, const int T, 
                      const float st_MX3[ NUM_NORMAL_STATES * (Q+1) * (T+1) ] );
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
void dp_matrix_Clear_X3 (const int Q, const int T,
                         float st_MX3[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                         float sp_MX[ NUM_SPECIAL_STATES * ((Q+1)+(T+1)) * 3 ],
                         int val);
void dp_matrix_Save (const int Q, const int T, 
                     const float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                     const float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ],
                     const char *_filename_);
void dp_matrix_trace_Save (const int Q, const int T, 
                        const float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                        const float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ],
                        const ALIGNMENT *aln,
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

#endif /* _TESTING_H */