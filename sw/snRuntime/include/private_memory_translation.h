#pragma once

#include "snrt.h"
#include "private_memory.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Translates a private memory address to a public memory address based on the cluster ID.
 *
 * This function takes a private memory address and a cluster ID, and translates the address
 * to a public memory address by adding an offset determined by the cluster ID.
 *
 * @param addr The private memory address to be translated.
 * @param cluster_id The cluster ID used to determine the offset for translation.
 * @return The translated public memory address.
 */

static inline uint32_t translate_to_public(uint32_t addr, uint32_t cluster_id) {
    return addr + (((uint32_t)cluster_id + 1)<<24);
}

#ifdef __cplusplus
}
#endif