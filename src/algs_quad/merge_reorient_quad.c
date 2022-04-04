/*******************************************************************************
 *  - FILE:      merge_reorient_quad.c
 *  - DESC:    Functions for merging multiple EDGEBOUND objects.
 *             Reorients from antidiagonal-wise to row-wise.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/* local imports */
#include "../objects/structs.h"
#include "../utilities/_utilities.h"
#include "../objects/_objects.h"

/* header */
#include "_algs_quad.h"
#include "merge_reorient_quad.h"

int EDGEBOUNDS_Merge_Reorient_Naive(const int Q,
                                    const int T,
                                    EDGEBOUNDS* edg_fwd,
                                    EDGEBOUNDS* edg_bck,
                                    EDGEBOUNDS* edg_diag,
                                    EDGEBOUNDS* edg_row,
                                    MATRIX_2D* cloud_MX) {
  /* merge edgebounds */
  MATRIX_2D_Fill(cloud_MX, 0);
  MATRIX_2D_Cloud_Fill(cloud_MX, edg_fwd, 1);
  MATRIX_2D_Cloud_Fill(cloud_MX, edg_bck, 2);

  /* reorient from diag-wise to row-wise */
  EDGEBOUNDS_Build_From_Cloud(Q, T, edg_row, cloud_MX, EDG_ROW);
  EDGEBOUNDS_Build_From_Cloud(Q, T, edg_diag, cloud_MX, EDG_DIAG);

  /* count number of cells covered by cloud */
  int num_cells = MATRIX_2D_Cloud_Count(cloud_MX);
  return num_cells;
}

int EDGEBOUNDS_Reorient_Naive(const int Q,
                              const int T,
                              EDGEBOUNDS* edg_dest,
                              EDGEBOUNDS* edg_src,
                              MATRIX_2D* cloud_MX) {
  int output_mode;
  if (edg_src->edg_mode == EDG_DIAG)
    output_mode = EDG_ROW;
  else
    output_mode = EDG_DIAG;

  MATRIX_2D_Fill(cloud_MX, 0);
  MATRIX_2D_Cloud_Fill(cloud_MX, edg_src, 1);
  EDGEBOUNDS_Build_From_Cloud(Q, T, edg_dest, cloud_MX, output_mode);
}

void EDGEBOUNDS_Build_From_Cloud(const int Q,
                                 const int T,
                                 EDGEBOUNDS* edg,
                                 MATRIX_2D* cloud_MX,
                                 int mode) {
  int d, i, j, k;
  int x, y1, y2;
  int d_st, d_end;
  int dim_min, dim_max;
  int le, re;
  int lb, rb;
  int num_cells;
  bool in_cloud;
  BOUND bnd;

  /* starts not in cloud */
  in_cloud = false;

  /* create edgebound in antidiagonal-wise order */
  if (mode == EDG_DIAG) {
    d_st = 0;
    d_end = Q + T;

    dim_min = MIN(Q, T);
    dim_max = MAX(Q, T);

    lb = 0;
    rb = 1;
    num_cells = 0;

    /* iterate through diags */
    for (d = d_st; d <= d_end; d++) {
      /* is dp matrix diagonal growing or shrinking? */
      if (d <= dim_min)
        num_cells++;
      if (d > dim_max)
        num_cells--;

      /* cells adjacent to previous row */
      lb = lb;
      rb = rb + 1;

      /* find diag cells that are inside matrix bounds */
      le = MAX(0, d - T);
      re = le + num_cells;

      /* for testing, bounds = edges */
      lb = MAX(lb, le);
      rb = MIN(rb, re);

      /* iterate through cells of diag */
      for (k = lb; k < rb; k++) {
        i = k;
        j = d - i;

        if (in_cloud) {
          /* end of bound, add bound to edgebound list */
          if (MX_2D(cloud_MX, i, j) <= 0) {
            in_cloud = false;
            y2 = k;
            BOUND bnd = (BOUND){x, y1, y2};
            EDGEBOUNDS_Pushback(edg, bnd);
          }
        } else {
          /* start of new bound */
          if (MX_2D(cloud_MX, i, j) > 0) {
            in_cloud = true;
            x = d;
            y1 = k;
          }
        }
      }

      /* if still in bound at end of diag, then at end of of bound, so add to edgebound list */
      if (in_cloud) {
        in_cloud = false;
        y2 = re;
        BOUND bnd = (BOUND){x, y1, y2};
        EDGEBOUNDS_Pushback(edg, bnd);
      }
    }
  }

  /* create edgebound in row-wise order */
  else if (mode == EDG_ROW) {
    /* iterate through rows */
    for (i = 0; i <= Q; i++) {
      /* iterate through cells in row */
      for (j = 0; j <= T; j++) {
        if (in_cloud) {
          /* end of bound, add bound to edgebound list */
          if (MX_2D(cloud_MX, i, j) <= 0) {
            in_cloud = false;
            y2 = j;

            BOUND bnd = (BOUND){x, y1, y2};
            EDGEBOUNDS_Pushback(edg, bnd);
          }
        } else {
          /* start of new bound */
          if (MX_2D(cloud_MX, i, j) > 0) {
            in_cloud = true;
            x = i;
            y1 = j;
            // printf("Starting [%d]:%d,(%d,?)...\n", EDGEBOUNDS_GetSize( edg ), x, y1);
          }
        }
      }

      /* if still in bound at end of row, then at end of of bound, so add to edgebound list */
      if (in_cloud) {
        in_cloud = false;
        y2 = T + 1;
        BOUND bnd = (BOUND){x, y1, y2};
        EDGEBOUNDS_Pushback(edg, bnd);
      }
    }
  }

  /* set edgebounds config */
  edg->Q = Q;
  edg->T = T;
  edg->edg_mode = mode;
}
