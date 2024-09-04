#pragma once

#include <stdint.h>

#include "snrt.h"


/**
 * @brief Initialize the event unit
 */
void teams_init(void);

/**
 * @brief send all workers in loop to exit()
 * @param global_core_idx global core index
 */
void teams_exit(uint32_t global_core_idx);

/**
 * @brief Enter the event unit loop, never exits
 *
 * @param global_core_idx global core index
 */
void teams_event_loop(uint32_t global_core_idx);

/**
 * @brief Set function to execute by `nthreads` number of threads
 * @details
 *
 * @param fn pointer to worker function to be executed
 * @param data pointer to function arguments
 * @param argc number of elements in data
 * @param nteams number of threads that have to execute this event
 */
int teams_dispatch_push(void (*fn)(void *, uint32_t), uint32_t argc, void *data,
                     uint32_t nteams);

/**
 * @brief wait for all workers to idle
 * @param global_core_idx global core index
 */
void teams_run_empty(uint32_t global_core_idx);

/**
 * @brief Acquires the event unit mutex, exits only on success
 */
void teams_mutex_lock();

/**
 * @brief Releases the acquired mutex
 */
void teams_mutex_release();

/**
 * Getters
 */
uint32_t teams_get_workers_in_loop();
uint32_t teams_get_workers_in_wfi();
