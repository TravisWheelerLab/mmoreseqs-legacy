/*******************************************************************************
 *  @file cloud_search.h
 *  @brief Testing for navigating through the matrices.
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/


#ifndef _TESTING_H
#define _TESTING_H

void fwd_test_cycle(int Q, int T,
                   float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                   float sp_MX[ NUM_NORMAL_STATES * (Q + 1) ],
                   TRACEBACK *tr);

void bck_test_cycle(int Q, int T,
                    float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                    float sp_MX[ NUM_NORMAL_STATES * (Q + 1) ],
                    TRACEBACK *tr);

int cloud_Fill(int Q, int T,
                float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                float sp_MX[ NUM_NORMAL_STATES * (Q + 1) ],
                EDGEBOUNDS* edg,
                float val, 
                int mode);

int cloud_Cell_Count(int Q, int T,
                   float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                   float sp_MX[ NUM_NORMAL_STATES * (Q + 1) ]);

#endif /* _TESTING_H */