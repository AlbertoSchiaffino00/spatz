#include "sync_teams.h"
#include <stdlib.h>

#include "printf.h"
#include "snrt.h"

inline void writeboh(uint32_t val, uintptr_t addr)
{
	asm volatile("sw %0, 0(%1)"
		     :
		     : "r"(val), "r"((volatile uint32_t *)addr)
		     : "memory");
}

static inline uint32_t readboh(const uintptr_t addr)
{
	uint32_t val;

	asm volatile("lw %0, 0(%1)"
		     : "=r"(val)
		     : "r"((const volatile uint32_t *)addr)
		     : "memory");
	return val;
}

typedef struct {
    uint32_t masters_in_loop;
    uint32_t exit_flag;
    uint32_t masters_mutex;
    uint32_t masters_wfi;
    struct {
        void (*fn)(void *, uint32_t);  // points to microtask wrapper
        void *data;
        uint32_t argc;
        uint32_t nteams;
        uint32_t fini_count;
    } e;
} teams_t;

volatile teams_t  teams_p_global = {0};

__thread volatile teams_t *teams_p;

//================================================================================
// prototypes
//================================================================================
static void wake_masters(void);
static void master_wfi(uint32_t global_core_idx);
static void wait_master_wfi(void);


//================================================================================
// public
//================================================================================
void teams_init(void) {
    if (snrt_global_core_idx() == 0) {
        // Allocate the eu struct in L1 for fast access
        // teams_p = snrt_l3alloc(sizeof(teams_t));
        // snrt_memset(&teams_p_global, 0, sizeof(teams_t));
        // store copy of teams_p on shared memory
        teams_p = &teams_p_global;

    } else if(snrt_cluster_core_idx() ==0 ){
        // while (!teams_p_global)
        //     ;
        teams_p = &teams_p_global;

    } // slaves of each cluster do not have access to the teams_p unit, for safety reasons
    snrt_cluster_hw_barrier();
}

void teams_exit(uint32_t global_core_idx) {

    // make sure queue is empty
    if (!teams_p->e.nteams) teams_run_empty(global_core_idx);
    // set exit flag and wake cores
    wait_master_wfi();

    teams_p->exit_flag = 1;

    wake_masters();
}

/**
 * @brief Return the number of workers currently present in the event loop
 */
uint32_t teams_get_masters_in_loop() {
    return __atomic_load_n(&teams_p->masters_in_loop, __ATOMIC_RELAXED);
}

/**
 * @brief Return the number of workers currently in the loop and waiting for
 * interrupt
 */
uint32_t teams_get_masters_in_wfi() {
    return __atomic_load_n(&teams_p->masters_wfi, __ATOMIC_RELAXED);
}



// Executed by all master threads, therefore (cluster_core_idx==0 && global_core_idx!=0) 

void teams_event_loop(uint32_t global_core_idx){
    
    uint32_t nteams;
    uint32_t cluster_idx = snrt_cluster_idx();

    // count number of workers in loop
    //__atomic_add_fetch(&teams_p->masters_in_loop, 1, __ATOMIC_RELAXED);
    teams_p->masters_in_loop += 1;

    // snrt_interrupt_enable(IRQ_M_SOFT);
    while (1) {


        // check for exit
        if (teams_p->exit_flag) {
            return;
        }


        if (cluster_idx < teams_p->e.nteams) {
            // make a local copy of nthreads to sync after work since the master
            // hart will reset eu_p->e.nthreads as soon as all workers finished
            // which might cause a race condition
            nteams = teams_p->e.nteams;

            // call    
            teams_p->e.fn(teams_p->e.data, teams_p->e.argc);

        }


        // enter wait for interrupt
        //__atomic_add_fetch(&teams_p->e.fini_count, 1, __ATOMIC_RELAXED);
        teams_p->e.fini_count += 1;
        master_wfi(global_core_idx);


    }
}

int teams_dispatch_push(void (*fn)(void *, uint32_t), uint32_t argc, void *data,
                     uint32_t nteams) {
    // wait for workers to be in wfi before manipulating the event struct
    wait_master_wfi();

    // fill queue
    teams_p->e.fn = fn;
    teams_p->e.data = data;
    teams_p->e.argc = argc;
    teams_p->e.nteams = nteams;


    return 0;
}

void teams_run_empty(uint32_t global_core_idx) {
    unsigned nfini, scratch;
    unsigned cluster_idx = snrt_cluster_idx();
    scratch = teams_p->e.nteams;


    if (!scratch) return;


    teams_p->e.fini_count = 0;
    if (scratch > 1) wake_masters();


    // Am i also part of the team?
    if (cluster_idx < teams_p->e.nteams) {
        // call  
        teams_p->e.fn(teams_p->e.data, teams_p->e.argc);

    }


    // wait for queue to be empty
    if (scratch > 1) {
        scratch = teams_get_masters_in_loop();
        while (__atomic_load_n(&teams_p->e.fini_count, __ATOMIC_RELAXED) !=
               scratch)
            ;
    }


    // stop workers from re-executing the task
    teams_p->e.nteams = 0;

}

/**
 * @brief Lock the event unit mutex
 */
inline void teams_mutex_lock() { snrt_mutex_lock(&teams_p->masters_mutex); }

/**
 * @brief Free the event unit mutex
 */
inline void teams_mutex_release() { snrt_mutex_release(&teams_p->masters_mutex); }


//================================================================================
// private
//================================================================================

static void wait_master_wfi(void) {

    uint32_t scratch = teams_p->masters_in_loop;
    while (__atomic_load_n(&teams_p->masters_wfi, __ATOMIC_RELAXED) != scratch)
        ;

}

static void wake_masters(void) {
    // wake all master cores except the main thread
    writeboh(1,0x40000804);
    writeboh(1,0x40000808);
    

}

static void master_wfi(uint32_t master_core_idx) {
    //__atomic_add_fetch(&teams_p->masters_wfi, 1, __ATOMIC_RELAXED);
    teams_p->masters_wfi += 1;

    snrt_wfi();

    // __atomic_add_fetch(&teams_p->masters_wfi, -1, __ATOMIC_RELAXED);
    teams_p->masters_wfi -= 1;

}

