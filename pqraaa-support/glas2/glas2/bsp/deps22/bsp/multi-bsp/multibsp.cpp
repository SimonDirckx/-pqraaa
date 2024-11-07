
#include "multibsp-internal.hpp"

#include <iostream>

namespace mcbsp {
	
	MultiBSP_program_algorithm_independent::MultiBSP_program_algorithm_independent( const size_t _p, const unsigned int _subtree_p, const unsigned int totalP ) : P( _p ), bsp_nleafs( totalP ), bsp_sleafs( _subtree_p ), bsp_lid( 0 ), bsp_P( 0 ), bsp_s( 0 ), state( UNDEFINED ) {}

	MultiBSP_program_algorithm_independent::MultiBSP_program_algorithm_independent() : MultiBSP_program_algorithm_independent( 0, 0, 0 ) { assert( false ); }

	MultiBSP_program_algorithm_independent::~MultiBSP_program_algorithm_independent() {}

	MultiBSP_template::~MultiBSP_template() {}

	void MultiBSP_program_algorithm_independent::monitor() {
		//spinlock until we can start initialisation
		while( state != UNDEFINED ) {}

		//initialise
		const int rc1 = pthread_mutex_init( &mutex, NULL );
		const int rc2 = pthread_cond_init( &condition, NULL );
		//DBG
		//std::cout << (&mutex) << " mutex initialised." << std::endl;
		if( rc1 != 0 || rc2 != 0 ) {
			bsp_abort(  "Error while initialising mutexes!\n" );
		}

		//prepare for idle state
		if( pthread_mutex_lock( &mutex ) != 0 ) {
			bsp_abort( "MultiBSP_program_algorithm_independent::monitor(): could not lock my own monitor!\n" );
		}

		//set suspended state
		state = SUSPENDED;

		//wake up when our state is set to RUNNING by our parent program
		while( state == SUSPENDED ) {
			//DBG
			//std::cout << "New program enters SUSPEND state." << std::endl;
			//otherwise, go to sleep
			pthread_cond_wait( &condition, &mutex );
		}

		//state should be RUNNING at this point
		assert( state == RUNNING );

		//DBG
		//std::cout << "New program enters RUNNING state." << std::endl;

		//cache nprocs/pids
		bsp_P = bsp_nprocs();
		bsp_s = bsp_pid();

		//go to user-defined function
		spmd();

		//automatic bsp_end
		state = EXITING;

		//let all siblings reach EXITING stage
		bsp_sync();

		//DBG
		//std::cout << "Thread is exiting" << std::endl;

		//wake up parent, allow them to reach exit state while we finish up
		wakeParent();

		//call recursive destruction
		finalize();

		//release lock
		if( pthread_mutex_unlock( &mutex ) != 0 ) {
			bsp_abort( "Could not release lock of my own monitor!\n" );
		}
		
		//let the callback function handle the cleaning of non-class member data;
		//only clean the MultiBSP-specific stuff here:
		const int rc7 = pthread_mutex_destroy( &mutex );
		const int rc8 = pthread_cond_destroy( &condition );
		if( rc7 != 0 || rc8 != 0 ) {
			std::cerr << "Non-fatal error: could not destroy monitor data structures upon MultiBSP program exit.\n";
		}

		//done
		return;
	}

	void MultiBSP_program_algorithm_independent::init() {
		//do leaf detection here
		if( P == 0 ) {
			//then no further BSP processes need be spawned
			return;
		}

		//proceed
		struct mcbsp_init_data * init = static_cast< struct mcbsp_init_data * >( malloc( sizeof(struct mcbsp_init_data) ) );

		//DBG
		//std::cerr << "Init struct created @ " << init << std::endl;

		pthread_attr_t attr;
#ifdef __MACH__
		const struct mcbsp_util_machine_info * const machine = mcbsp_internal_getMachineInfo();
#elif defined _WIN32
		DWORD_PTR mask;
#else
		cpu_set_t mask;
#endif
		//initialise init object
		init->spmd = NULL;
		init->bsp_program = NULL;
		init->argc  = 0;
		init->argv  = NULL;
		init->threads = static_cast< pthread_t* >( malloc( P*sizeof(pthread_t) ) );
		init->P     = P;
		init->abort = false;
		init->ended = false;
		init->sync_entry_counter = 0;
		init->sync_exit_counter  = 0;
		pthread_mutex_init( &(init->mutex), NULL );
#ifdef MCBSP_USE_SPINLOCK
		init->condition     = static_cast< unsigned char * >(malloc( P * sizeof( unsigned char ) ));
		init->mid_condition = static_cast< unsigned char * >(malloc( P * sizeof( unsigned char ) ));
		if( init->condition == NULL || init->mid_condition == NULL ) {
			std::cerr << "Could not allocate spinlock conditions!\n";
			mcbsp_util_fatal();
		}
		for( size_t s = 0; s < ((size_t)P); ++s ) {
			init->condition[ s ] = init->mid_condition[ s ] = 0;
		}
#else
		pthread_cond_init ( &(init->condition), NULL );
		pthread_cond_init ( &(init->mid_condition), NULL );
#endif
		mcbsp_util_address_table_initialise( &(init->global2local), ((size_t)P) );
		init->threadData = static_cast< struct mcbsp_thread_data * * >(malloc( ((size_t)P) * sizeof( struct mcbsp_thread_data * ) ));
		init->prev_data  = static_cast< struct mcbsp_thread_data * >(pthread_getspecific( mcbsp_internal_thread_data ));
		init->tagSize = 0;

		//DBG
		//std::cout << "Requesting pinning for " << P << " threads." << std::endl;

		//derive pinning
		struct mcbsp_util_pinning_info pinning_info = mcbsp_util_pinning( ((size_t)P), mcbsp_internal_getMachineInfo(),
			init->prev_data == NULL ?
			NULL :
			&(init->prev_data->machine_partition) );
		size_t *pinning = pinning_info.pinning;

		//sanity checks
		if( pinning == NULL ) {
			std::cerr <<  "Could not get a valid pinning!\n";
			mcbsp_util_fatal();
		}

#ifdef MCBSP_SHOW_PINNING
		//report pinning:
		std::cout << "Info: pinning used is";
		for( size_t s=0; s < ((size_t)P); ++s ) {
			std::cout << " " << ((unsigned long int)pinning[ s ]);
		}
		std::cout << "\n";
#endif

		//spawn P threads.
		for( size_t s = ((size_t)P) - 1; s < ((size_t)P); --s ) {
			//get new thread_data object, with initialised plain old data
			struct mcbsp_thread_data *thread_data = mcbsp_internal_allocate_thread_data( init, s );
			//provide the submachine of this thread (note this is also plain old data in thread_data)
			if( pinning_info.partition.top == 0 ) { //manual partitioning; no submachines available. Derive a simple one.
				thread_data->machine_partition.start = thread_data->machine_partition.end = pinning_info.pinning[ s ];
				thread_data->machine_partition.stride = 1;
			} else { //get submachines from the pinning info
				const struct mcbsp_util_machine_partition * const partitions_array = (struct mcbsp_util_machine_partition*)(pinning_info.partition.array);
				thread_data->machine_partition = partitions_array[ s ];
			}

			//spawn new threads for each sub-process
			//create POSIX threads attributes
			//(currently for pinning, if supported)
			pthread_attr_init( &attr );
#ifdef _WIN32
			mask = (DWORD_PTR)1;
			mask <<= pinning[ s ];
#elif !defined(__MACH__)
			CPU_ZERO( &mask );
			CPU_SET ( pinning[ s ], &mask );
			pthread_attr_setaffinity_np( &attr, sizeof( cpu_set_t ), &mask );
#endif
			//get corresponding instance
			MultiBSP_program_algorithm_independent * const target_instance = getChildPointer( s );
			//save thread-local data
			target_instance->thread_data = static_cast< void * >(thread_data);
			//spawn the actual thread
			const int pthr_rval = pthread_create( &(init->threads[ s ]), &attr, &MultiBSP_program_algorithm_independent::callback, getChildPointer( s ) );
			if( pthr_rval != 0 ) {
				std::cerr << "Could not spawn new thread (" << strerror( pthr_rval ) << ")!\n";
				mcbsp_util_fatal();
			}
			//store reference to child with ID 0
			if( s == 0 ) {
				nested_process = init->threads[ 0 ];
			}
#ifdef _WIN32
			//do after-creation pinning
			SetThreadAffinityMask( pthread_getw32threadhandle_np( init->threads[ s ] ), mask );
#elif defined __MACH__
			//set after-creation affinities
			thread_port_t osx_thread = pthread_mach_thread_np( init->threads[ s ] );
			struct thread_affinity_policy ap;
			switch( machine->affinity ) {
				case SCATTER:
				{
					//Affinity API release notes do not specify whether 0 is a valid tag, or in fact equal to NULL; so 1-based to be sure
					ap.affinity_tag = s + 1;
					break;
				}
				case COMPACT:
				{
					ap.affinity_tag = 1;
					break;
				}
				case MANUAL:
				{
					ap.affinity_tag = machine->manual_affinity[ s ];
					break;
				}
				default:
				{
					std::cerr << "Unhandled affinity type for Mac OS X!\n";
					mcbsp_util_fatal();
				}
			}
			thread_policy_set( osx_thread, THREAD_AFFINITY_POLICY, (thread_policy_t)&ap, THREAD_AFFINITY_POLICY_COUNT );
#endif
			//destroy attributes object
			pthread_attr_destroy( &attr );
		}

		//free pinning_info
		free( pinning );
		mcbsp_util_stack_destroy( &(pinning_info.partition) );
	}

	void * MultiBSP_program_algorithm_independent::callback( void * const p_instance ) {
		//get instance pointer
		MultiBSP_program_algorithm_independent &instance = * static_cast< MultiBSP_program_algorithm_independent * >( p_instance );
		//get thread-local data
		struct mcbsp_thread_data * data = static_cast< struct mcbsp_thread_data * >( instance.thread_data );

		//reallocate thread-local data and initialise communication queues
		instance.thread_data = (data = mcbsp_internal_initialise_thread_data( data ));

		//provide a link back from the initialising process' data
		data->init->threadData[ data->bsp_id ] = data;

		//store thread-local data
		if( pthread_setspecific( mcbsp_internal_thread_data, data ) != 0 ) {
			std::cerr << "Could not store thread local data!\n";
			mcbsp_util_fatal();
		}
#ifdef __MACH__
		//get rights for accessing Mach's timers
		const kern_return_t rc1 = host_get_clock_service( mach_host_self(), SYSTEM_CLOCK, &(data->clock) );
		if( rc1 != KERN_SUCCESS ) {
			std::cerr << "Could not access the Mach system timer (" << mach_error_string( rc1 ) << ")\n";
			mcbsp_util_fatal();
		}

		//record start time
		const kern_return_t rc2 = clock_get_time( data->clock, &(data->start) );
		if( rc2 != KERN_SUCCESS ) {
			std::cerr << "Could not get start time (" << mach_error_string( rc2 ) << ")\n";
			mcbsp_util_fatal();
		}
#elif defined _WIN32
		//record start time
		QueryPerformanceCounter( &(data->start) );
		//record frequency
		QueryPerformanceFrequency( &(data->frequency) );
#else
		//record start time
		clock_gettime( CLOCK_MONOTONIC, &(data->start) );
#endif
		//initialise recursively to spawn child threads also
		instance.init();
		//continue with the instance monitor daemon code
		instance.monitor();
		//store thread-local data
		if( pthread_setspecific( mcbsp_internal_thread_data, NULL ) != 0 ) {
			std::cerr << "Could not store thread local data!\n";
			mcbsp_util_fatal();
		}
		//cache init data
		struct mcbsp_init_data * const init = data->init;
		//cache PID
		const size_t pid = data->bsp_id;
		//destroy thread data
		mcbsp_internal_destroy_thread_data( data );
		//if we are PID 0, do a synchronised exit
		if( pid == 0 ) {
			//first wait for other threads to exit
			for( size_t s = 1; s < init->P; ++s ) {
				if( pthread_join( init->threads[ s ], NULL ) != 0 ) {
					std::cerr << "Error: MultiBSP_program_algorithm_independent::callback call to pthread_join failed!" << std::endl;
					exit( EXIT_FAILURE );
				}
			}
			//it's on us to clear the init struct
			pthread_mutex_destroy( &(init->mutex) );
#ifdef MCBSP_USE_SPINLOCK
			free( init->condition );
			free( init->mid_condition );
#else
			pthread_cond_destroy( &(init->    condition) );
			pthread_cond_destroy( &(init->mid_condition) );
#endif
			mcbsp_util_address_table_destroy( &(init->global2local) );
			free( init->threadData );
			free( init->threads );
			//destroy init data
			free( init );
		}
		//exit cleanly
		return NULL;
	}

};

