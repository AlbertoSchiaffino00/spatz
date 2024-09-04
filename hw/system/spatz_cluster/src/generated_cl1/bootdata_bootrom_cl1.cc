// Copyright 2021 ETH Zurich and University of Bologna.
// Solderpad Hardware License, Version 0.51, see LICENSE for details.
// SPDX-License-Identifier: SHL-0.51

#include <stdint.h>

// The boot data generated along with the system RTL.
struct BootData {
    uint32_t boot_addr;
    uint32_t core_count;
    uint32_t hartid_base;
    uint32_t tcdm_start;
    uint32_t tcdm_size;
    uint32_t tcdm_offset;
    uint64_t global_mem_start;
    uint64_t global_mem_end;
};

extern "C" const BootData BOOTDATA = {.boot_addr = 0x1000,
                           .core_count = 2,
                           .hartid_base = 16,
                           .tcdm_start = 0x51000000,
                           .tcdm_size = 0x20000,
                           .tcdm_offset = 0x0,
                           .global_mem_start = 0x78000000,
                           .global_mem_end = 0x78400000};