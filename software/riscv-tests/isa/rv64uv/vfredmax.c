// Copyright 2022 ETH Zurich and University of Bologna.
// Solderpad Hardware License, Version 0.51, see LICENSE for details.
// SPDX-License-Identifier: SHL-0.51
//
// Author: Xiaorui Yin <yinx@student.ethz.ch>
// Date: 2022/05/03

#include "float_macros.h"
#include "vector_macros.h"

// Naive test
void TEST_CASE1(void) {
  VSET(16, e16, m8);
  // 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8
  VLOAD_16(v16, 0x3c00, 0x4000, 0x4200, 0x4400, 0x4500, 0x4600, 0x4700, 0x4800,
           0x3c00, 0x4000, 0x4200, 0x4400, 0x4500, 0x4600, 0x4700, 0x4800);
  VLOAD_16(v24, 0x3c00);
  asm volatile("vfredmax.vs v8, v16, v24");
  VCMP_U16(1, v8, 0x4800);

  VSET(16, e32, m8);
  VLOAD_32(v16, 0x3F800000, 0x40000000, 0x40400000, 0x40800000, 0x40A00000,
           0x40C00000, 0x40E00000, 0x41000000, 0x3F800000, 0x40000000,
           0x40400000, 0x40800000, 0x40A00000, 0x40C00000, 0x40E00000,
           0x41000000);
  VLOAD_32(v24, 0x3F800000);
  asm volatile("vfredmax.vs v8, v16, v24");
  VCMP_U32(2, v8, 0x41000000);

#if ELEN == 64
  VSET(16, e64, m8);
  VLOAD_64(v16, 0x3FF0000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
  VLOAD_64(v24, 0x3FF0000000000000);
  asm volatile("vfredmax.vs v8, v16, v24");
  VCMP_U64(3, v8, 0x4020000000000000);
#endif

  // Super long vector length
  VSET(64, e32, m8);
  VLOAD_32(
      v16, 0x3F800000, 0x40000000, 0x40400000, 0x40800000, 0x40A00000,
      0x40C00000, 0x40E00000, 0x41000000, 0x3F800000, 0x40000000, 0x40400000,
      0x40800000, 0x40A00000, 0x40C00000, 0x40E00000, 0x41000000, 0x3F800000,
      0x40000000, 0x40400000, 0x40800000, 0x40A00000, 0x40C00000, 0x40E00000,
      0x41000000, 0x3F800000, 0x40000000, 0x40400000, 0x40800000, 0x40A00000,
      0x40C00000, 0x40E00000, 0x41000000, 0x3F800000, 0x40000000, 0x40400000,
      0x40800000, 0x40A00000, 0x40C00000, 0x40E00000, 0x41000000, 0x3F800000,
      0x40000000, 0x40400000, 0x40800000, 0x40A00000, 0x40C00000, 0x40E00000,
      0x41000000, 0x3F800000, 0x40000000, 0x40400000, 0x40800000, 0x40A00000,
      0x40C00000, 0x40E00000, 0x41000000, 0x3F800000, 0x40000000, 0x40400000,
      0x40800000, 0x40A00000, 0x40C00000, 0x40E00000, 0x41000000);
  VLOAD_32(v24, 0x3F800000);
  asm volatile("vfredmax.vs v8, v16, v24");
  VCMP_U32(4, v8, 0x41000000);
}

// Masked naive test
void TEST_CASE2(void) {
  VSET(16, e16, m8);
  VLOAD_8(v0, 0xaa, 0x55);
  // 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8
  VLOAD_16(v16, 0x3c00, 0x4000, 0x4200, 0x4400, 0x4500, 0x4600, 0x4700, 0x4800,
           0x3c00, 0x4000, 0x4200, 0x4400, 0x4500, 0x4600, 0x4700, 0x4800);
  VLOAD_16(v24, 0x3c00);
  asm volatile("vfredmax.vs v8, v16, v24, v0.t");
  VCMP_U16(5, v8, 0x4800);

  VSET(16, e32, m8);
  VLOAD_8(v0, 0xaa, 0x55);
  VLOAD_32(v16, 0x3F800000, 0x40000000, 0x40400000, 0x40800000, 0x40A00000,
           0x40C00000, 0x40E00000, 0x41000000, 0x3F800000, 0x40000000,
           0x40400000, 0x40800000, 0x40A00000, 0x40C00000, 0x40E00000,
           0x41000000);
  VLOAD_32(v24, 0x3F800000);
  asm volatile("vfredmax.vs v8, v16, v24, v0.t");
  VCMP_U32(6, v8, 0x41000000);

#if ELEN == 64
  VSET(16, e64, m8);
  VLOAD_8(v0, 0xaa, 0x55);
  VLOAD_64(v16, 0x3FF0000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
  VLOAD_64(v24, 0x3FF0000000000000);
  asm volatile("vfredmax.vs v8, v16, v24, v0.t");
  VCMP_U64(7, v8, 0x4020000000000000);
#endif
}

// Are we respecting the undisturbed tail policy?
void TEST_CASE3(void) {
  VSET(16, e16, m8);
  // 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8
  VLOAD_16(v16, 0x3c00, 0x4000, 0x4200, 0x4400, 0x4500, 0x4600, 0x4700, 0x4800,
           0x3c00, 0x4000, 0x4200, 0x4400, 0x4500, 0x4600, 0x4700, 0x4800);
  VLOAD_16(v24, 0x3c00, 0x4000, 0x4200, 0x4400, 0x4500, 0x4600, 0x4700, 0x4800,
           0x3c00, 0x4000, 0x4200, 0x4400, 0x4500, 0x4600, 0x4700, 0x4800);
  VLOAD_16(v8, 0x3c00, 0x4000, 0x4200, 0x4400, 0x4500, 0x4600, 0x4700, 0x4800,
           0x3c00, 0x4000, 0x4200, 0x4400, 0x4500, 0x4600, 0x4700, 0x4800);
  asm volatile("vfredmax.vs v8, v16, v24");
  VCMP_U16(8, v8, 0x4800, 0x4000, 0x4200, 0x4400, 0x4500, 0x4600, 0x4700,
           0x4800, 0x3c00, 0x4000, 0x4200, 0x4400, 0x4500, 0x4600, 0x4700,
           0x4800);

  VSET(16, e32, m8);
  VLOAD_32(v16, 0x3F800000, 0x40000000, 0x40400000, 0x40800000, 0x40A00000,
           0x40C00000, 0x40E00000, 0x41000000, 0x3F800000, 0x40000000,
           0x40400000, 0x40800000, 0x40A00000, 0x40C00000, 0x40E00000,
           0x41000000);
  VLOAD_32(v24, 0x3F800000, 0x40000000, 0x40400000, 0x40800000, 0x40A00000,
           0x40C00000, 0x40E00000, 0x41000000, 0x3F800000, 0x40000000,
           0x40400000, 0x40800000, 0x40A00000, 0x40C00000, 0x40E00000,
           0x41000000);
  VLOAD_32(v8, 0x3F800000, 0x40000000, 0x40400000, 0x40800000, 0x40A00000,
           0x40C00000, 0x40E00000, 0x41000000, 0x3F800000, 0x40000000,
           0x40400000, 0x40800000, 0x40A00000, 0x40C00000, 0x40E00000,
           0x41000000);
  asm volatile("vfredmax.vs v8, v16, v24");
  VCMP_U32(9, v8, 0x41000000, 0x40000000, 0x40400000, 0x40800000, 0x40A00000,
           0x40C00000, 0x40E00000, 0x41000000, 0x3F800000, 0x40000000,
           0x40400000, 0x40800000, 0x40A00000, 0x40C00000, 0x40E00000,
           0x41000000);

#if ELEN == 64
  VSET(16, e64, m8);
  VLOAD_64(v16, 0x3FF0000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
  VLOAD_64(v24, 0x3FF0000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
  VLOAD_64(v8, 0x3FF0000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
  asm volatile("vfredmax.vs v8, v16, v24");
  VCMP_U64(10, v8, 0x4020000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
#endif
}

// Odd number of elements, undisturbed policy
void TEST_CASE4(void) {
#if ELEN == 64
  VSET(1, e64, m8);
  VLOAD_64(v16, 0x3FF0000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
  VLOAD_64(v24, 0x3FF0000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
  VLOAD_64(v8, 0x3FF0000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
  asm volatile("vfredmax.vs v8, v16, v24");
  VCMP_U64(11, v8, 0x3FF0000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);

  VSET(3, e64, m8);
  VLOAD_64(v16, 0x3FF0000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
  VLOAD_64(v24, 0x3FF0000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
  VLOAD_64(v8, 0x3FF0000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
  asm volatile("vfredmax.vs v8, v16, v24");
  VCMP_U64(12, v8, 0x4008000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3ff0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);

  VSET(7, e64, m8);
  VLOAD_64(v16, 0x3FF0000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
  VLOAD_64(v24, 0x3FF0000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
  VLOAD_64(v8, 0x3FF0000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
  asm volatile("vfredmax.vs v8, v16, v24");
  VCMP_U64(13, v8, 0x401C000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3ff0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);

  VSET(15, e64, m8);
  VLOAD_64(v16, 0x3FF0000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
  VLOAD_64(v24, 0x3FF0000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
  VLOAD_64(v8, 0x3FF0000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
  asm volatile("vfredmax.vs v8, v16, v24");
  VCMP_U64(14, v8, 0x4020000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
#endif
}

// Odd number of elements, undisturbed policy, and mask
void TEST_CASE5(void) {
  VSET(7, e16, m8);
  VLOAD_8(v0, 0x00, 0xff);
  // 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8
  VLOAD_16(v16, 0x3c00, 0x4000, 0x4200, 0x4400, 0x4500, 0x4600, 0x4700, 0x4800,
           0x3c00, 0x4000, 0x4200, 0x4400, 0x4500, 0x4600, 0x4700, 0x4800);
  VLOAD_16(v24, 0x3c00, 0x4000, 0x4200, 0x4400, 0x4500, 0x4600, 0x4700, 0x4800,
           0x3c00, 0x4000, 0x4200, 0x4400, 0x4500, 0x4600, 0x4700, 0x4800);
  VLOAD_16(v8, 0x3c00, 0x4000, 0x4200, 0x4400, 0x4500, 0x4600, 0x4700, 0x4800,
           0x3c00, 0x4000, 0x4200, 0x4400, 0x4500, 0x4600, 0x4700, 0x4800);
  asm volatile("vfredmax.vs v8, v16, v24, v0.t");
  VCMP_U16(15, v8, 0x3c00, 0x4000, 0x4200, 0x4400, 0x4500, 0x4600, 0x4700,
           0x4800, 0x3c00, 0x4000, 0x4200, 0x4400, 0x4500, 0x4600, 0x4700,
           0x4800);

  VSET(1, e32, m8);
  VLOAD_8(v0, 0xff, 0x00);
  VLOAD_32(v16, 0x3F800000, 0x40000000, 0x40400000, 0x40800000, 0x40A00000,
           0x40C00000, 0x40E00000, 0x41000000, 0x3F800000, 0x40000000,
           0x40400000, 0x40800000, 0x40A00000, 0x40C00000, 0x40E00000,
           0x41000000);
  VLOAD_32(v24, 0x3F800000, 0x40000000, 0x40400000, 0x40800000, 0x40A00000,
           0x40C00000, 0x40E00000, 0x41000000, 0x3F800000, 0x40000000,
           0x40400000, 0x40800000, 0x40A00000, 0x40C00000, 0x40E00000,
           0x41000000);
  VLOAD_32(v8, 0x3F800000, 0x40000000, 0x40400000, 0x40800000, 0x40A00000,
           0x40C00000, 0x40E00000, 0x41000000, 0x3F800000, 0x40000000,
           0x40400000, 0x40800000, 0x40A00000, 0x40C00000, 0x40E00000,
           0x41000000);
  asm volatile("vfredmax.vs v8, v16, v24, v0.t");
  VCMP_U32(16, v8, 0x3F800000, 0x40000000, 0x40400000, 0x40800000, 0x40A00000,
           0x40C00000, 0x40E00000, 0x41000000, 0x3F800000, 0x40000000,
           0x40400000, 0x40800000, 0x40A00000, 0x40C00000, 0x40E00000,
           0x41000000);

#if ELEN == 64
  VSET(3, e64, m8);
  VLOAD_8(v0, 0xaa, 0x55);
  VLOAD_64(v16, 0x3FF0000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
  VLOAD_64(v24, 0x3FF0000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
  VLOAD_64(v8, 0x3FF0000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
  asm volatile("vfredmax.vs v8, v16, v24, v0.t");
  VCMP_U64(17, v8, 0x4000000000000000, 0x4000000000000000, 0x4008000000000000,
           0x4010000000000000, 0x4014000000000000, 0x4018000000000000,
           0x401C000000000000, 0x4020000000000000, 0x3FF0000000000000,
           0x4000000000000000, 0x4008000000000000, 0x4010000000000000,
           0x4014000000000000, 0x4018000000000000, 0x401C000000000000,
           0x4020000000000000);
#endif
}

int main(void) {
  INIT_CHECK();
  enable_vec();
  enable_fp();

  TEST_CASE1();
  // TEST_CASE2();
  TEST_CASE3();
  TEST_CASE4();
  // TEST_CASE5();

  EXIT_CHECK();
}
