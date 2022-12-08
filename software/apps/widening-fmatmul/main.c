// Copyright 2021 ETH Zurich and University of Bologna.
//
// SPDX-License-Identifier: Apache-2.0
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Author: Domenic Wüthrich, ETH Zurich

#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include "data/data_gemm.h"
#include "kernel/widening-fmatmul.c"
#include "printf.h"
#ifdef MEMPOOL
#include "alloc.h"
#include "runtime.h"
#include "synchronization.h"
#endif

// Define Matrix dimensions:
// C = AB with A=[MxN], B=[NxP], C=[MxP]
#ifndef MATRIX_DIM
#define MATRIX_DIM 256
#endif
#ifndef KERNEL_M
#define KERNEL_M 2
#endif
#ifndef KERNEL_P
#define KERNEL_P 2
#endif

#define M MATRIX_DIM
#define N MATRIX_DIM
#define P MATRIX_DIM

__fp16 *a;
__fp16 *b;
__fp16 *c;

// Initialize the matrices
void init_matrix(double *dst, const double *src, const unsigned int len) {
  for (unsigned int i = 0; i < len; ++i) {
    dst[i] = src[i];
  }
}

// Verify the matrices
int verify_matrix(__fp16 *matrix, const __fp16 *checksum,
                  const unsigned int num_rows, const unsigned int num_columns) {
  for (unsigned int i = 0; i < num_rows; ++i) {
    float sum = 0;
    for (unsigned int j = 0; j < num_columns; ++j) {
      sum += (float)matrix[i * num_columns + j];
    }

    float diff = sum - (float)checksum[i];
    if (diff < 0)
      diff = -diff;
    if (diff > 0.001) {
      return i == 0 ? -1 : (int)i;
    }
  }
  return 0;
}

void print_matrix(__fp16 const *matrix, unsigned int num_rows,
                  unsigned int num_columns) {
  printf("0x%8X\n", (unsigned int)matrix);
  for (unsigned int i = 0; i < num_rows; ++i) {
    for (unsigned int j = 0; j < num_columns; ++j) {
      printf("%5u ", (unsigned int)matrix[i * num_columns + j]);
    }
    printf("\n");
  }
}

int main() {
  const unsigned int num_cores = mempool_get_core_count();
  const unsigned int cores_per_group = num_cores / NUM_GROUPS;
  const unsigned int cid = mempool_get_core_id();
  const unsigned int core_gid = cid % cores_per_group;
  const unsigned int gid = cid / cores_per_group;

  const unsigned int active_groups = 1;
  const unsigned int active_cores = cores_per_group * active_groups;
  const unsigned int is_core_active = cid < active_cores;

  const unsigned int measure_iterations = 1;

  unsigned int timer_start, timer_end, timer;
  unsigned int row_start, row_end;

  unsigned int m_start, m_end;
  unsigned int p_start, p_end;
  unsigned int vl;
  unsigned int dim;
  unsigned int kernel_size;

  // Initialize MemPool
  mempool_init(cid, num_cores);

  // Initialize multicore barrier
  mempool_barrier_init(cid);

  // Allocate the matrices in the local tile
  if (cid == 0) {
    a = (__fp16 *)domain_malloc(get_alloc_tile(0), M * N * sizeof(__fp16));
    b = (__fp16 *)domain_malloc(get_alloc_tile(0), N * P * sizeof(__fp16));
    c = (__fp16 *)domain_malloc(get_alloc_tile(0), M * P * sizeof(__fp16));
  }

  // Reset timer
  timer = (unsigned int)-1;

  // Set matrix dimension
  dim = MATRIX_DIM;
  kernel_size = KERNEL_M;
  vl = KERNEL_P;

  // Can every core execute its desired kernel?
  if ((dim * dim) / (kernel_size * vl) < active_cores)
    return -1;
  // Does the vl fit inside the dim
  if (vl > dim)
    return -2;

  // Block dimension of group
  const unsigned int dim_group = dim / active_groups;
  // Number of parallel cores in m direction
  const unsigned int split_m_count = dim_group / kernel_size;

  if (split_m_count < cores_per_group) {
    // Split P dimension up
    const unsigned int split_p_count = cores_per_group / split_m_count;
    p_start = dim / split_p_count * (core_gid % split_p_count);
    p_end = dim / split_p_count * ((core_gid % split_p_count) + 1);
    m_start = dim_group * gid + kernel_size * (core_gid / split_p_count);
    m_end = dim_group * gid + kernel_size * (core_gid / split_p_count + 1);
  } else {
    // Work over complete P dimension
    p_start = 0;
    p_end = dim;
    m_start = dim_group * gid + (dim_group / cores_per_group) * core_gid;
    m_end = dim_group * gid + (dim_group / cores_per_group) * (core_gid + 1);
  }

  // Wait for all cores to finish
  mempool_barrier(num_cores);

  // Initialize matrices
  const unsigned int cores_per_row = active_cores / dim;
  if (dim < active_cores) {
    row_start = cid / cores_per_row;
    row_end = cid / cores_per_row + 1;
  } else {
    row_start = dim / num_cores * cid;
    row_end = dim / num_cores * (cid + 1);
  }

  if (cid == 0) {
    init_matrix((double*)a, (const double*)gemm_A_dram, dim * dim * sizeof(__fp16)/ sizeof(double));
    init_matrix((double*)b, (const double*)gemm_B_dram, dim * dim * sizeof(__fp16)/ sizeof(double));
    init_matrix((double*)c, (const double*)gemm_C_dram, dim * dim * sizeof(__fp16)/ sizeof(double));
  }

  // Wait for all cores to finish
  mempool_barrier(num_cores);

  for (unsigned int i = 0; i < measure_iterations; ++i) {
    // Calculate matmul
    if (is_core_active) {
      // Start timer
      timer_start = mempool_get_timer();

      if (kernel_size == 2) {
        matmul_2xVL(c, a, b, m_start, m_end, dim, dim, p_start, p_end, vl);
      } else if (kernel_size == 4) {
        matmul_4xVL(c, a, b, m_start, m_end, dim, dim, p_start, p_end, vl);
      } else if (kernel_size == 8) {
        matmul_8xVL(c, a, b, m_start, m_end, dim, dim, p_start, p_end, vl);
      } else {
        return -2;
      }
    }

    // Wait for all cores to finish matmul
    mempool_barrier(num_cores);

    // End timer and check if new best runtime
    timer_end = mempool_get_timer();
    unsigned int timer_temp = timer_end - timer_start;
    if (cid == 0) {
      if (timer_temp < timer) {
        timer = timer_temp;
      }
    }
  }

  // Check and display results
  if (cid == 0) {
    unsigned int performance = 1000 * 2 * dim * dim * dim / timer;
    unsigned int utilization = performance / (2 * active_cores * 2 * N_IPU);

    printf("\n----- (%dx%d) sdotp fmatmul -----\n", dim, dim);
    printf("The execution took %u cycles.\n", timer);
    printf("The performance is %u OP/1000cycle (%u%%o utilization).\n",
           performance, utilization);
  }

  if (cid == 0) {
    int error = verify_matrix(c, (const __fp16 *)gemm_checksum, dim, dim);

    if (error != 0) {
      printf("Error core %d: c[%d]=%u\n", cid, error, (int)c[error]);
      return error;
    }
  }

  // Free the matrices
  if (cid == 0) {
    domain_free(get_alloc_tile(0), a);
    domain_free(get_alloc_tile(0), b);
    domain_free(get_alloc_tile(0), c);
  }

  // Wait for core 0 to finish displaying results
  mempool_barrier(num_cores);

  return 0;
}
