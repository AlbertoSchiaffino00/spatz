// Copyright 2020 ETH Zurich and University of Bologna.
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0

#include "stdint.h"

static inline void writew(uint32_t val, uintptr_t addr)
{
	asm volatile("sw %0, 0(%1)"
		     :
		     : "r"(val), "r"((volatile uint32_t *)addr)
		     : "memory");
}


static inline uint32_t readw(const uintptr_t addr)
{
	uint32_t val;

	asm volatile("lw %0, 0(%1)"
		     : "=r"(val)
		     : "r"((const volatile uint32_t *)addr)
		     : "memory");
	return val;
}

int main() { 

    //should return the starting code address
    return readw(0x50020058);
}
