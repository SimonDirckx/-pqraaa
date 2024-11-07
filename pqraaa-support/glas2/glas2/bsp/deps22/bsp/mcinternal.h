/*
 * Copyright (c) 2012 by Albert-Jan N. Yzelman
 *
 * This file is part of MulticoreBSP in C --
 *        a port of the original Java-based MulticoreBSP.
 *
 * MulticoreBSP for C is distributed as part of the original
 * MulticoreBSP and is free software: you can redistribute
 * it and/or modify it under the terms of the GNU Lesser 
 * General Public License as published by the Free Software 
 * Foundation, either version 3 of the License, or 
 * (at your option) any later version.
 * MulticoreBSP is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the 
 * implied warranty of MERCHANTABILITY or FITNESS FOR A 
 * PARTICULAR PURPOSE. See the GNU Lesser General Public 
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General 
 * Public License along with MulticoreBSP. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef _H_MCINTERNAL
#define _H_MCINTERNAL

#include "mcbsp-internal.h"
#include "mcutil.h"

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <assert.h>
#include <pthread.h>
#include <stdbool.h>

#ifdef __MACH__
 #include <mach/mach.h>
 #include <mach/clock.h>
 #include <mach/mach_error.h>
#endif

#ifdef _WIN32
 #include <windows.h>
#endif

/**
 * Initialisation struct.
 * 
 * Contains information necessary to start an SPMD
 * section, and contains fields necessary for global
 * actions during the corresponding SPMD section.
 *
 * Can be stored at two locations:
 * 1: at the initialising thread of an SPMD section,
 *    when the SPMD section has not yet been started;
 * 2: as a pointer through thread-local data once the
 *    SPMD section has started.
 */
struct mcbsp_init_data {

	/** User-definied SPMD entry-point. */
	void (*spmd)(void);

	/** 
	 * In case of a call from the C++ wrapper,
	 * pointer to the user-defined BSP_program
	 */
	void *bsp_program;

	/** Passed argc from bsp_init. */
	int argc;

	/** Passed argv from bsp_init. */
	char **argv;

	/** Threads corresponding to this BSP program. */
	pthread_t *threads;

	/** Number of processors involved in this run. */
	size_t P;

	/** Whether the BSP program should be aborted. */
	volatile bool abort;

	/** Whether the BSP program has ended. */
	bool ended;

	/**
	 * Barrier counter for this BSP execution.
	 * Synchronises synchronisation entry.
	 */
	size_t sync_entry_counter;

	/**
	 * Barrier counter for this BSP program.
	 * Synchronises synchronisation exit.
	 */
	size_t sync_exit_counter;

	/** Mutex used for critical sections (such as synchronisation). */
	pthread_mutex_t mutex;

#ifdef MCBSP_USE_SPINLOCK
	/** Condition used for critical sections. */
	unsigned char * condition;

	/** Condition used for critical sections. */
	unsigned char * mid_condition;
#else
	/** Condition used for critical sections. */
	pthread_cond_t condition;

	/** Condition used for critical sections. */
	pthread_cond_t mid_condition;
#endif

	/** Address table used for inter-thread communication. */
	struct mcbsp_util_address_table global2local;

	/** Pointers to all thread-local data, as needed for communication. */
	struct mcbsp_thread_data **threadData;

	/** Stores any previous thread-local data. Used for nested runs. */
	struct mcbsp_thread_data *prev_data;

	/** Currently active tag size. */
	size_t tagSize;

};

/**
 * Thread-local data.
 *
 * Contains data local to each thread in an SPMD section.
 */
struct mcbsp_thread_data {

	/** Thread-local BSP id. */
	size_t bsp_id;

	/**
	 * Stores the tag size to become active after the
	 * next synchronisation.
	 */
	size_t newTagSize;

	/**
	 * Counts the maximum number of registered variables
	 * at any one time.
	 */
	size_t localC;

	/** Initialisation data. */
	struct mcbsp_init_data * init;

	/**
	 * The communication queues used for bsp_get
	 * requests.
	 */
	struct mcbsp_util_stack * request_queues;

	/**
	 * The communication queues used for all DRMA
	 * and BSMP communication.
	 */
	struct mcbsp_util_stack * queues;

	/** 
	 * The communication queues used for bsp_hpsend
	 * requests. Kept separately from bsp_hpget and
	 * bsp_hpput queues since they are handled quite
	 * differently.
	 */
	struct mcbsp_util_stack * hpsend_queues;

	/** Local address to global variable map. */
	struct mcbsp_util_address_map local2global;

	/**
	 * The communication queue for hpget and hpput
	 * requests. There is only a single queue per
	 * thread, since only the active threads knows
	 * when all requests to-be are made and thus
	 * when the communication can start.
	 * This is what makes hp-primitives unsafe; 
	 * overlapping or neighbouring memory areas can
	 * cause data races and false sharing conditions.
	 * The main advantage of using the hp queue is
	 * that no payload data will be cached.
	 */
	struct mcbsp_util_stack hpdrma_queue;

	/**
	 * Keeps track of which global variables are removed
	 * (as per bsp_pop_reg).
	 */
	struct mcbsp_util_stack removedGlobals;

	/**
	 * Keeps track which globals will be removed before
	 * the next superstep arrives.
	 */
	struct mcbsp_util_stack localsToRemove;

	/** The BSMP incoming message queue. */
	struct mcbsp_util_stack bsmp;

	/** The push request queue. */
	struct mcbsp_util_stack localsToPush;

	/**
	 * Caches the mapSearch results for use with
	 * bsp_push_reg requests.
	 */
	struct mcbsp_util_stack globalsToPush;

	/** Stores the machine part allocated to this thread. Used for nested runs. */
	struct mcbsp_util_machine_partition machine_partition;

#ifdef __MACH__
	/** OS port for getting timings */
	clock_serv_t clock;

	/** 
	 * Stores the start-time of this thread.
	 * (Mach-specific (OS X) timespec)
	 */
	mach_timespec_t start;
#elif defined _WIN32
	/** 
	 * Stores the start-time of this thread.
	 * (Windows-specific timespec)
	 */
	LARGE_INTEGER start;

	/**
	 * Frequency of the Windows high-
	 * resultion timer.
	 */
	LARGE_INTEGER frequency;
#else
	/** Stores the start-time of this thread. */
	struct timespec start;
#endif
};

/**
 *  A DRMA communication request for `get'-requests.
 *  @see bsp_get
 */
struct mcbsp_get_request {

	/** Source */
	void * source;

	/** Destination */
	void * destination;

	/** Length */
	size_t length;

};

/**
 *  A generic BSP communication message.
 *  Handles both `put' and regular BSMP requests.
 *  @see bsp_put
 *  @see bsp_send
 *  @see mcbsp_get_request (used for `get'-requests)
 *  @see mcbsp_hpsend_request (used for high-performance BSMP) 
 */
struct mcbsp_message {
	
	/** Destination */
	void * destination;

	/** Length */
	size_t length;

};

/**
 *  A DRMA communication request for hpput
 *  and hpget requests.
 */
struct mcbsp_hp_request {

	/** Source */
	const void * source;

	/** Destination */
	void * destination;

	/** Length */
	size_t length;

};

/**
 *  A high-performance (non-buffering) BSMP message request.
 *  @see bsp_hpsend.
 */
struct mcbsp_hpsend_request {

	/** Payload source */
	const void * payload;

	/** Tag source */
	const void * tag;

	/** Payload length */
	size_t payload_size;

};

/** Struct corresponding to a single push request. */
struct mcbsp_push_request {

	/** The local address to push. */
	void * address;

	/** The memory range to register. */
	size_t size;

};

/** Per-initialising thread initialisation data. */
extern pthread_key_t mcbsp_internal_init_data;

/** Per-thread data. */
extern pthread_key_t mcbsp_internal_thread_data;

/** Per-thread machine info. */
extern pthread_key_t mcbsp_local_machine_info;
/**
 * Per-thread machine info.
 * Necessary for nested MulticoreBSP SPMD calls.
 */
extern pthread_key_t mcbsp_internal_init_data;

/** Whether mcbsp_internal_init_data is initialised. */
extern bool mcbsp_internal_keys_allocated;

/** 
 * Contorls thread-safe singleton access to
 * mcbsp_internal_init_data.
 * Required for safe initialisation.
 */
extern pthread_mutex_t mcbsp_internal_keys_mutex;

/**
 * Performs a BSP intialisation using the init struct supplied.
 * The construction of this struct differs when called from C code,
 * or when called from the C++ wrapper.
 */
void bsp_init_internal( struct mcbsp_init_data * const initialisationData );

/**
 * Checks if everything is all right to start an SPMD program.
 *
 * If there is an error, execution is stopped using mcbsp_util_fatal.
 *
 * @return NULL if the call to bsp_begin was valid, but no further
 * 		action is required; or a pointer to an initialisation
 * 		struct when the SPMD program is yet to be spawned.
 */
struct mcbsp_init_data * bsp_begin_check( void );

/**
 * Gets the machine info currently active for this
 * MulticoreBSP session.
 *
 * @return A pointer to the current machine info struct.
 */
struct mcbsp_util_machine_info *mcbsp_internal_getMachineInfo( void );

/**
 * Singleton thread-safe allocator for mcbsp_internal_init_data.
 */
void mcbsp_internal_check_keys_allocated( void );

/**
 * Entry-point of MulticoreBSP threads. Initialises internals
 * and then executes the user-defined SPMD program.
 *
 * @param p The passed POSIX threads init data.
 */
void* mcbsp_internal_spmd( void *p );

/**
 * Checks if a abort has been requested, and if so,
 * exits the current thread.
 */
void mcbsp_internal_check_aborted( void );

#ifdef MCBSP_USE_SPINLOCK
/**
 * Alias for mcbsp_internal_syncWithCondition using
 * the standard init.condition
 *
 * @param init   Pointer to the BSP init corresponding
 *               to our current SPMD group.
 * @param bsp_id The unique ID number corresponding to
 *               the thread that calls this sync function.
 */
void mcbsp_internal_sync( struct mcbsp_init_data * const init, const size_t bsp_id );
#else
/**
 * Alias for mcbsp_internal_syncWithCondition using
 * the standard init.condition
 *
 * @param init Pointer to the BSP init corresponding
 *             to our current SPMD group.
 */
void mcbsp_internal_sync( struct mcbsp_init_data * const init );
#endif

#ifdef MCBSP_USE_SPINLOCK
/**
 * Implements an actual synchronisation. This is
 * used within bsp_sync, but also in bsp_begin to
 * ensure all is initialised before getting into
 * a real communication-enabled sync (necessary
 * to support high-performant put, get, or send 
 * operations).
 *
 * @param init Pointer to the BSP init corresponding
 *             to our current SPMD group.
 * @param condition The condition to use with
 *                  this particular sync.
 * @param bsp_id The unique ID number corresponding 
 *               to the thread that calls this sync
 *               function.
 */
void mcbsp_internal_syncWithCondition( struct mcbsp_init_data * const init,
	volatile unsigned char * const condition, const size_t bsp_id );
#else
/**
 * Implements an actual synchronisation. This is
 * used within bsp_sync, but also in bsp_begin to
 * ensure all is initialised before getting into
 * a real communication-enabled sync (necessary
 * to support high-performant put, get, or send 
 * operations).
 *
 * @param init Pointer to the BSP init corresponding
 *             to our current SPMD group.
 * @param condition The condition to use with
 *                  this particular sync.
 */
void mcbsp_internal_syncWithCondition( struct mcbsp_init_data * const init,
	pthread_cond_t * const condition );
#endif

/**
 * Common part executed by all BSP primitives when in SPMD part.
 * This version assumes local thread data used by BSP remains
 * unchanged.
 */ 
const struct mcbsp_thread_data * mcbsp_internal_const_prefunction( void );

/**
 * Common part executed by all BSP primitives when in SPMD part.
 * This is the non-const version of mcbsp_internal_const_prefunction.
 */
struct mcbsp_thread_data * mcbsp_internal_prefunction( void );

/** Default SPMD function to call. */
int main();

/**
 * Allcates and initialises the `plain old data' parts of the
 * mcbsp_thread_data struct, using the supplied parameters.
 * @param init Pointer to the initialising process' init data.
 * @param s    Process ID of the thread corresponding to this data.
 * @return The newly allocated thread data.
 */
struct mcbsp_thread_data * mcbsp_internal_allocate_thread_data( struct mcbsp_init_data * const init, const size_t s );

/**
 * Finishes initialisation of a mcbsp_thread_data struct. The
 * plain-old-data fields of a given thread_data are copied,
 * and the remaining data fields are properly initialised.
 * The old thread_data struct will be freed.
 *
 * @param data The initial data to copy.
 * @return New process-local thread data.
 */
struct mcbsp_thread_data * mcbsp_internal_initialise_thread_data( struct mcbsp_thread_data * const data );

/**
 * Deallocates a fully initialised mcbsp_thread_data struct.
 *
 * @param data The thread data to destroy.
 */
void mcbsp_internal_destroy_thread_data( struct mcbsp_thread_data * const data );

#endif

