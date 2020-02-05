/*******************************************************************************
 *  @file edgebounds_obj.h
 *  @brief EDGEBOUNDS Object functions.
 *
 *  @author Dave Rich
 *  @bug Lots.
 *******************************************************************************/


#ifndef _EDGEBOUNDS_OBJ_H
#define _EDGEBOUNDS_OBJ_H

EDGEBOUNDS *edgebounds_Create();

void edgebounds_Init(EDGEBOUNDS **edg);

void edgebounds_Destroy(EDGEBOUNDS *edg);

void edgebounds_Add(EDGEBOUNDS *edg,
                    BOUND bnd);

void edgebounds_Resize(EDGEBOUNDS *edg);

void edgebounds_Print(EDGEBOUNDS *edg);
void bound_Print(BOUND bnd);

void edgebounds_Save(EDGEBOUNDS *edg,
                     const char *_filename_);

void edgbounds_Sort(EDGEBOUNDS *edg);

void edgebounds_Merge(int Q, int T,
                     EDGEBOUNDS *edg_fwd,
                     EDGEBOUNDS *edg_bck,
                     EDGEBOUNDS *edg_new);

void edgebounds_Reorient(int Q, int T,
                         EDGEBOUNDS *edg_diag,
                         EDGEBOUNDS *edg_row);

int edgebounds_Merge_Reorient_Cloud(EDGEBOUNDS*edg_fwd,
                                    EDGEBOUNDS*edg_bck,
                                    EDGEBOUNDS*edg_new,
                                    int Q, int T,
                                    float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                                    float sp_MX[ NUM_NORMAL_STATES * (Q + 1) ]);

void edgebounds_Reorient_Cloud( EDGEBOUNDS*edg_old,
                               EDGEBOUNDS*edg_new,
                               int Q, int T,
                               float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                               float sp_MX[ NUM_NORMAL_STATES * (Q + 1) ], 
                               int old_mode, int new_mode);

void edgebounds_Merge_Cloud(EDGEBOUNDS*edg_1,
                            EDGEBOUNDS*edg_2,
                            EDGEBOUNDS*edg_res,
                            int Q, int T,
                            float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                            float sp_MX[ NUM_NORMAL_STATES * (Q + 1) ],
                            int mode);

void edgebounds_Build_From_Cloud( EDGEBOUNDS*edg,
                                 int Q, int T,
                                 float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                                 int mode);

#endif /* _EDGEBOUNDS_OBJ_H */