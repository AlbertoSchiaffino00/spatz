// Copyright 2021 ETH Zurich and University of Bologna.
// Solderpad Hardware License, Version 0.51, see LICENSE for details.
// SPDX-License-Identifier: SHL-0.51

#include <tb_lib.hh>

namespace sim {

const BootData BOOTDATA = {.boot_addr = 0x1000,
                           .core_count = 2,
                           .hartid_base = 18,
                           .tcdm_start = 0x52000000,
                           .tcdm_size = 0x20000,
                           .tcdm_offset = 0x0,
                           .global_mem_start = 0x78000000,
                           .global_mem_end = 0x78400000};

}  // namespace sim
