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

//for documentation see the .h file:
#include "mcutil.h"

size_t mcbsp_util_detect_hardware_threads( void ) {
#ifdef _WIN32
	SYSTEM_INFO system_info;
	GetSystemInfo( &system_info );
	return (size_t) system_info.dwNumberOfProcessors;
#else
	return (size_t) sysconf( _SC_NPROCESSORS_ONLN );
#endif
}

struct mcbsp_util_machine_info *mcbsp_util_createMachineInfo( void ) {
	struct mcbsp_util_machine_info * MCBSP_MACHINE_INFO = malloc( sizeof( struct mcbsp_util_machine_info ) );
	if( MCBSP_MACHINE_INFO == NULL ) {
		fprintf( stderr, "Error: could not allocate new machine info struct!\n" );
		mcbsp_util_fatal();
	}
	//set defaults
	MCBSP_MACHINE_INFO->Tset   = false;
	MCBSP_MACHINE_INFO->Cset   = false;
	MCBSP_MACHINE_INFO->Aset   = false;
	MCBSP_MACHINE_INFO->TpCset = false;
	MCBSP_MACHINE_INFO->TNset  = false;
	MCBSP_MACHINE_INFO->Pset   = false;
	MCBSP_MACHINE_INFO->Rset   = false;
	MCBSP_MACHINE_INFO->UTset  = false;
	MCBSP_MACHINE_INFO->threads 		= 0;
	MCBSP_MACHINE_INFO->affinity 		= MCBSP_DEFAULT_AFFINITY;
	MCBSP_MACHINE_INFO->cores 		= 0;
	MCBSP_MACHINE_INFO->threads_per_core 	= MCBSP_DEFAULT_THREADS_PER_CORE;
	MCBSP_MACHINE_INFO->thread_numbering	= MCBSP_DEFAULT_THREAD_NUMBERING;
	MCBSP_MACHINE_INFO->manual_affinity	= NULL;
	MCBSP_MACHINE_INFO->num_reserved_cores  = 0;
	MCBSP_MACHINE_INFO->reserved_cores	= NULL;
	MCBSP_MACHINE_INFO->unused_threads_per_core = 0;

	//keep track whether we have set machine settings
	//that don't have straightforward defaults
	bool Tset = false;
	bool Cset = false;
	//check for machine.info file to parse and get settings from
	FILE *fp = fopen ( "machine.info", "r" ) ;
	if( fp != NULL ) {
		//file exists, allocate buffers
		char LINE_BUFFER[ 255 ];
		char RECURSIVE_LINE_BUFFER[ 255 ];
		char key[ 255 ];
		//read line-by-line
		while( fgets( LINE_BUFFER, 255, fp ) != NULL ) {
			//get key, and process if value was not already changed from code
			if( sscanf( LINE_BUFFER, "%s %[^\n]", key, RECURSIVE_LINE_BUFFER ) == 2 ) {
				//check for cores key
				if( !MCBSP_MACHINE_INFO->Cset && strcmp( key, "cores" ) == 0 ) {
					unsigned long int value = 0;
					if( sscanf( LINE_BUFFER, "%s %lu", key, &value ) != 2 ) {
						fprintf( stderr, "Warning: parsing error while processing `cores' key in machine.info at: %s", LINE_BUFFER );
						assert( Cset == false );
					} else {
						MCBSP_MACHINE_INFO->cores = (size_t) value;
						Cset = true;
					}
				//check for threads key
				} else if( !MCBSP_MACHINE_INFO->Tset && strcmp( key, "threads" ) == 0 ) {
					unsigned long int value = 0;
					if( sscanf( LINE_BUFFER, "%s %lu", key, &value ) != 2 ) {
						fprintf( stderr, "Warning: parsing error while processing `threads' key in machine.info at: %s", LINE_BUFFER );
						assert( Tset == false );
					} else {
						MCBSP_MACHINE_INFO->threads = (size_t) value;
						Tset = true;
					}
				//check for threads_per_core key
				} else if( !MCBSP_MACHINE_INFO->TpCset && strcmp( key, "threads_per_core" ) == 0 ) {
					unsigned long int value = 0;
					bool success = true;
					if( sscanf( LINE_BUFFER, "%s %lu", key, &value ) != 2 ) {
						fprintf( stderr, "Warning: parsing error while processing `threads_per_core' key (reverting to default) in machine.info at: %s",
							LINE_BUFFER );
						success = false;
					} else {
						MCBSP_MACHINE_INFO->threads_per_core = value;
					}
					if( !success ) {
						MCBSP_MACHINE_INFO->threads_per_core = MCBSP_DEFAULT_THREADS_PER_CORE;
					}
				//check for unused_threads_per_core key
				} else if( !MCBSP_MACHINE_INFO->UTset && strcmp( key, "unused_threads_per_core" ) == 0 ) {
					unsigned long int value = 0;
					bool success = true;
					if( sscanf( LINE_BUFFER, "%s %lu", key, &value ) != 2 ) {
						fprintf( stderr, "Warning: parsing error while processing `unused_threads_per_core' key in machine.info at: %s", LINE_BUFFER );
						success = false;
					} else {
						MCBSP_MACHINE_INFO->unused_threads_per_core = value;
					}
					if( !success ) {
						MCBSP_MACHINE_INFO->unused_threads_per_core = 0;
					}	
				//check for threads_numbering key
				} else if( !MCBSP_MACHINE_INFO->TNset && strcmp( key, "thread_numbering" ) == 0 ) {
					char value[ 255 ];
					bool success = true;
					if( sscanf( LINE_BUFFER, "%s %s", key, value ) != 2 ) {
						fprintf( stderr, "Warning: parsing error while processing `thread_numbering' key (reverting to default) in machine.info at: %s", LINE_BUFFER );
						success = false;
					} else {
						if( strcmp( value, "consecutive" ) == 0 ) {
							MCBSP_MACHINE_INFO->thread_numbering = CONSECUTIVE;
						} else if( strcmp( value, "wrapped" ) == 0 ) {
							MCBSP_MACHINE_INFO->thread_numbering = WRAPPED;
						} else {
							fprintf( stderr, "Warning: unkown value for the `thread_numbering' key in machine.info (%s); reverting to default value.\n", value );
							success = false;
						}
					}
					if( !success ) {
						MCBSP_MACHINE_INFO->thread_numbering = MCBSP_DEFAULT_THREAD_NUMBERING;
					}
				//check for affinity key
				} else if( !MCBSP_MACHINE_INFO->Aset && strcmp( key, "affinity" ) == 0 ) {
					char value[ 255 ];
					bool success = true;
					if( sscanf( LINE_BUFFER, "%s %s", key, value ) != 2 ) {
						fprintf( stderr, "Warning: parsing error while processing `affinity' key in machine.info (reverting to default) at: %s", LINE_BUFFER );
						success = false;
					} else {
						if( strcmp( value, "scatter" ) == 0 ) {
							MCBSP_MACHINE_INFO->affinity = SCATTER;
						} else if( strcmp( value, "compact" ) == 0 ) {
							MCBSP_MACHINE_INFO->affinity = COMPACT;
						} else if( strcmp( value, "manual" ) == 0 ) {
							MCBSP_MACHINE_INFO->affinity = MANUAL;
						} else {
							fprintf( stderr, "Warning: unkown value for the `affinity' key in machine.info (%s); reverting to default value.\n", value );
							success = false;
						}
					}
					//if an error occurred, revert to default
					if( !success ) {
						MCBSP_MACHINE_INFO->affinity = MCBSP_DEFAULT_AFFINITY;
					}
				//echeck for pinning key
				} else if( !MCBSP_MACHINE_INFO->Pset && strcmp( key, "pinning" ) == 0 ) {
					//if P was not set yet, set P to the maximum and assume user will supply complete pinning list
					if( !Tset ) {
						MCBSP_MACHINE_INFO->threads = mcbsp_util_detect_hardware_threads();
						Tset = true;
					}
					//allocate manual pinning array
					MCBSP_MACHINE_INFO->manual_affinity = malloc( (MCBSP_MACHINE_INFO->threads) * sizeof( size_t ) );
					//prepare to read pinning list of size P, and assume success
					size_t cur = 0;
					bool success = true;
					//prepare reading of vector
					unsigned long int value;
					int number_read;
					//recursive read until input line is exhausted
					do {
						number_read = sscanf( RECURSIVE_LINE_BUFFER, "%lu %[^\n]", &value, RECURSIVE_LINE_BUFFER );
						if( number_read < 1 ) {
							fprintf( stderr, "Warning: error in parsing a value for `pinning' keyword in machine.info at: %s\n", RECURSIVE_LINE_BUFFER );
							success = false;
						} else {
							MCBSP_MACHINE_INFO->manual_affinity[ cur++ ] = (size_t)value;
						}
					} while( success && number_read == 2 );
					//check whether the input line had exactly enough values
					if( cur < MCBSP_MACHINE_INFO->threads ) {
						fprintf( stderr, "Warning: not enough values for `pinning' keyword (read %lu) in machine.info at: %s", 
							(unsigned long int)cur, LINE_BUFFER );
						success = false;
					} else if( cur > MCBSP_MACHINE_INFO->threads ) {
						fprintf( stderr, "Warning: too many values for `pinning' keyword in machine.info (remainder will be ignored) at: %s", LINE_BUFFER );
					}
					//if reading of the pinning list failed, fall back to sane defaults
					if( !success ) {
						fprintf( stderr, "Warning: setting manual pinning failed, reverting to default affinity strategy if necessary.\n" );
						free( MCBSP_MACHINE_INFO->manual_affinity );
						MCBSP_MACHINE_INFO->manual_affinity = NULL;
						if( MCBSP_MACHINE_INFO->affinity == MANUAL ) {
							MCBSP_MACHINE_INFO->affinity = MCBSP_DEFAULT_AFFINITY;
						}
					}
				} else if( !MCBSP_MACHINE_INFO->Rset && strcmp( key, "reserved_cores" ) == 0 ) {
					//first parse size
					MCBSP_MACHINE_INFO->num_reserved_cores = 0;
					unsigned long int value;
					int number_read;
					bool success = true;
					do { //(same strategy as for pinning, see above)
						number_read = sscanf( RECURSIVE_LINE_BUFFER, "%lu %[^\n]", &value, RECURSIVE_LINE_BUFFER );
						++(MCBSP_MACHINE_INFO->num_reserved_cores);
					} while( number_read == 2 );
					if( MCBSP_MACHINE_INFO->num_reserved_cores == 0 ) {
						fprintf( stderr, "Warning: number of reserved cores (%lu) does not make sense.\n",
							(unsigned long int)(MCBSP_MACHINE_INFO->num_reserved_cores)
						 );
						success = false;
					}
					if( success ) {
						//allocate the list
						MCBSP_MACHINE_INFO->reserved_cores = malloc( (MCBSP_MACHINE_INFO->num_reserved_cores) * sizeof( size_t ) );
						if( MCBSP_MACHINE_INFO->reserved_cores == NULL ) {
							fprintf( stderr, "Warning: could not allocate reserved_cores array; reverting back to defaults (no reserved cores).\n" );
							success = false;
						}
					}
					if( success ) {
						//parse the list
						sscanf( LINE_BUFFER, "%s %[^\n]", key, RECURSIVE_LINE_BUFFER ); //reset RECURSIVE_LINE_BUFFER
						size_t cur = 0;
						do {
							number_read = sscanf( RECURSIVE_LINE_BUFFER, "%lu %[^\n]", &value, RECURSIVE_LINE_BUFFER );
							if( number_read < 1 ) {
								fprintf( stderr, "Warning: error in parsing a value for `reserved_cores' keyword in machine.info at: %s\n", RECURSIVE_LINE_BUFFER );
								success = false; 
							} else {
								MCBSP_MACHINE_INFO->reserved_cores[ cur++ ] = (size_t)value;
							}
						} while( success && number_read == 2 );
					}
					if( !success ) {
						if( MCBSP_MACHINE_INFO->reserved_cores != NULL ) {
							free( MCBSP_MACHINE_INFO->reserved_cores );
						}
						MCBSP_MACHINE_INFO->num_reserved_cores = 0;
						MCBSP_MACHINE_INFO->reserved_cores = NULL;
					}
				} else {
					fprintf( stderr, "Warning: unknown key in machine.info (%s)\n", key );
				}
			} else {
				fprintf( stderr, "Warning: error while reading machine.info at: %s", LINE_BUFFER );
			}
		} //end while-not-(EOF-or-I/O-error) loop
		//close file handle
		if( fclose( fp ) != 0 ) {
			fprintf( stderr, "Warning: error while closing the machine.info file." );
		}
	} //end read file branch
	const size_t syscores = mcbsp_util_detect_hardware_threads();
	//check #cores was set, if not, supply default (maximum available CPUs on this system / threads_per_core)
	if( !Cset ) {
		MCBSP_MACHINE_INFO->cores = syscores / MCBSP_MACHINE_INFO->threads_per_core + (syscores % MCBSP_MACHINE_INFO->threads_per_core > 0 ? 1 : 0);
	}
	//check if P (#threads) was set, if not, supply default (maximum available CPUs on this system)
	if( !Tset ) {
		MCBSP_MACHINE_INFO->threads = (MCBSP_MACHINE_INFO->cores - MCBSP_MACHINE_INFO->num_reserved_cores) *
						(MCBSP_MACHINE_INFO->threads_per_core - MCBSP_MACHINE_INFO->unused_threads_per_core);
	}
	//user did not set cores nor threads explicitly (yet)
	//(previous usage was internal)
	Cset = false;
	Tset = false;

	//sanity checks, but no hard fails
	if( MCBSP_MACHINE_INFO->cores <= MCBSP_MACHINE_INFO->num_reserved_cores ) {
		fprintf( stderr,
			"Warning: number of reserved cores (%ld) is larger than or equal to the number of available cores (%lu). Attempting to continue anyway...\n",
			(unsigned long int)(MCBSP_MACHINE_INFO->num_reserved_cores),
			(unsigned long int)(MCBSP_MACHINE_INFO->cores)
		);
	}
	if( MCBSP_MACHINE_INFO->threads_per_core <= MCBSP_MACHINE_INFO->unused_threads_per_core ) {
		fprintf( stderr,
			"Warning: number of unused threads per core (%ld) leaves no available threads (since threads_per_core equals %ld). Attempting to continue anyway...\n",
			(unsigned long int)(MCBSP_MACHINE_INFO->unused_threads_per_core),
			(unsigned long int)(MCBSP_MACHINE_INFO->threads_per_core)
		);
	}
	if( (MCBSP_MACHINE_INFO->cores - MCBSP_MACHINE_INFO->num_reserved_cores) * 
		(MCBSP_MACHINE_INFO->threads_per_core - MCBSP_MACHINE_INFO->unused_threads_per_core) < MCBSP_MACHINE_INFO->threads ) {
		fprintf( stderr,
			"Warning: number of threads (%lu) is set higher than the number of available cores (%lu-%lu) times the number of available threads per core (%lu-%lu).",
			(unsigned long int)(MCBSP_MACHINE_INFO->threads),
			(unsigned long int)(MCBSP_MACHINE_INFO->cores),
			(unsigned long int)(MCBSP_MACHINE_INFO->num_reserved_cores),
			(unsigned long int)(MCBSP_MACHINE_INFO->threads_per_core),
			(unsigned long int)(MCBSP_MACHINE_INFO->unused_threads_per_core)
		);
		fprintf( stderr, " Attempting to continue anyway...\n" );
	}
	//done
	return MCBSP_MACHINE_INFO;
}

void mcbsp_util_destroyMachineInfo( void * machine_info ) {
	//if no NULL
	if( machine_info != NULL ) {
		//cast pointer to actual data
		struct mcbsp_util_machine_info * MCBSP_MACHINE_INFO =
			(struct mcbsp_util_machine_info *) machine_info;
		//if a manual affinity is defined
		if( MCBSP_MACHINE_INFO->manual_affinity != NULL ) {
			//delete it
			free( MCBSP_MACHINE_INFO->manual_affinity );
			MCBSP_MACHINE_INFO->manual_affinity = NULL;
		}
		//if reserved cores/threads were defined
		if( MCBSP_MACHINE_INFO->reserved_cores != NULL ) {
			free( MCBSP_MACHINE_INFO->reserved_cores );
			MCBSP_MACHINE_INFO->num_reserved_cores = 0;
			MCBSP_MACHINE_INFO->reserved_cores = NULL;
		}
		//delete it
		free( MCBSP_MACHINE_INFO );
		MCBSP_MACHINE_INFO = NULL;
	}
}

int mcbsp_util_integer_compare( const void * a, const void * b ) {
	//normal comparison, plus guard against overflow
	const size_t x = *((size_t*)a);
	const size_t y = *((size_t*)b);
	const signed char d1 = x > y;
	const signed char d2 = x < y;
	return d1 - d2;
}

size_t mcbsp_util_log2( size_t x ) {
	size_t ret = 0; 
	while( x >>= 1 ) //while there are 1-bits after a right shift,
		++ret;   //increase the log by one
	return ret;
}

size_t mcbsp_util_sort_unique_integers( size_t * const array, const size_t length, const size_t lower_bound, const size_t upper_bound ) {
	if( length == 0 ) return 0; //an empty array is already sorted
	//sanity checks
	if( array == NULL ) {
		fprintf( stderr, "Error: NULL array passed to mcbsp_util_sort_unique_integers!\n" );
		mcbsp_util_fatal();
	}
	if( lower_bound >= upper_bound ) {
		fprintf( stderr, "Error: invalid value range passed to mcbsp_util_sort_unique_integers (%lu,%lu)!\n",
			(unsigned long int)lower_bound, (unsigned long int)upper_bound );
		mcbsp_util_fatal();
	}
	//allocate helper-array
	bool * counting_sort = malloc( (upper_bound - lower_bound ) * sizeof( bool ) );
	if( counting_sort == NULL ) {
		fprintf( stderr, "Error: could not allocate helper-array for counting sort!\n" );
		mcbsp_util_fatal();
	}
	//initialise helper-array
	for( size_t s = 0; s < upper_bound - lower_bound; ++s )
		counting_sort[ s ] = false; //mark integer s as unused
	//do sort, phase I
	bool warned = false;
	for( size_t s = 0; s < length; ++s ) {
		const size_t value = array[ s ];
		if( value < lower_bound || value >= upper_bound ) {
			fprintf( stderr, "Error: array to sort contained out-of-range element (mcbsp_util_sort_unique_integers in mcutil.c)!\n" );
			mcbsp_util_fatal();
		} else {
			const size_t toUpdate = value - lower_bound;
			if( !warned && counting_sort[ toUpdate ] ) {
				fprintf( stderr, "Warning: array to sort contains duplicate values (mcbsp_util_sort_unique_integers in mcutil.c)!\n" );
				warned = true;
			} else {
				counting_sort[ array[ s ] - lower_bound ] = true; //mark the value of array[ s ] as used
			}
		}
	}
	//do sort, phase II
	size_t cur = 0;
	for( size_t s = 0; s < upper_bound - lower_bound; ++s ) {
		if( counting_sort[ s ] ) //then the value s+lower_bound was taken
			array[ cur++ ] = s + lower_bound;
	}
	//free helper array
	free( counting_sort );
	//done
	return cur;
}

bool mcbsp_util_contains( const size_t * const array, const size_t value, const size_t lo, const size_t hi ) {
	//sanity
	if( lo >= hi ) {
		fprintf( stderr, "Error: bounds given to mcbsp_util_contains (mcutil.c) do not make sense (%lu,%lu)!\n",
			(unsigned long int)lo, (unsigned long int)hi );
		mcbsp_util_fatal();
	}

	//check lower bound
	if( array[ lo ] == value )
		return true;

	//check upper bound
	if( array[ hi - 1 ] == value )
		return true;

	//get middle value
	const size_t mid = ( lo + hi ) / 2;

	//check for end of recursion
	if( mid == lo )
		return false;
	
	//recurse
	if( array[ mid ] > value )
		return mcbsp_util_contains( array, value, lo, mid );
	else
		return mcbsp_util_contains( array, value, mid, hi );

}

struct mcbsp_util_pinning_info mcbsp_util_pinning( const size_t P, struct mcbsp_util_machine_info * const machine, const struct mcbsp_util_machine_partition * const super_machine ) {
	//allocate return variable
	struct mcbsp_util_pinning_info ret;

	//initialise return pointer
	ret.pinning = NULL;
	
	//check boundary condition
	if( P == 0 ) {
		mcbsp_util_stack_initialise( &(ret.partition), sizeof( struct mcbsp_util_machine_partition ) );
		return ret;
	}

	//check machine info
	mcbsp_util_check_machine_info( machine );

	//available number of cores for this MulticoreBSP run
	const size_t available_cores   = machine->cores - machine->num_reserved_cores;
	const size_t available_threads = machine->Tset ? machine->threads : available_cores * (machine->threads_per_core-machine->unused_threads_per_core);

	//check if there are sufficient threads
	if( P > available_threads ) {
		fprintf( stderr, "Warning: %lu threads requested but only %lu available. If this was on purpose, please supply a proper manual affinity and increase the maximum threads usable by MulticoreBSP manually (see the documentation for details)\n", (unsigned long int)P, (unsigned long int)(available_threads) );
	}

	//handle the case of manual affinity
	if( machine->affinity == MANUAL ) {
		//We do create local copies so to not conflict with possible nested bsp_begin()'s.
		//(Not for a bit of extra data locality: it is read only once.)
		ret.pinning = malloc( P * sizeof( size_t ) );
		if( ret.pinning == NULL ) {
			fprintf( stderr, "Error: could not allocate local pinning array!\n" );
			mcbsp_util_fatal();
		}
		if( machine->manual_affinity == NULL ) {
			fprintf( stderr, "No manual affinity array provided, yet a manual affinity strategy was set.\n" );
			mcbsp_util_fatal();
		}
		for( size_t s = 0; s < P; ++s ) {
			ret.pinning[ s ] = machine->manual_affinity[ s ];
		}
		mcbsp_util_stack_initialise( &(ret.partition), sizeof( struct mcbsp_util_machine_partition ) );
	} else {
		//handle systematic affinities; first derive bottom-level machine descriptors
		if( super_machine == NULL ) {
			ret.partition = mcbsp_util_partition( mcbsp_util_derive_partition( machine ), P, machine->affinity );
		} else {
			ret.partition = mcbsp_util_partition( *super_machine, P, machine->affinity );
		}

		//sanity checks
		if( mcbsp_util_stack_empty( &(ret.partition) ) ) {
			fprintf( stderr, "Error: did not get a valid machine partitioning during call to mcbsp_util_pinning.\n" );
			mcbsp_util_fatal();
		}
		if( ret.partition.top != P ) {
			fprintf( stderr, "Error: expected %ld parts but mcbsp_util_partition returned only %ld.\n", P, ret.partition.top );
			mcbsp_util_fatal();
		}
		//get standard pinning
		const struct mcbsp_util_machine_partition * const partitions_array = (struct mcbsp_util_machine_partition*)(ret.partition.array);
		ret.pinning = mcbsp_util_derive_standard_pinning( partitions_array, P );
		//filter pinning for reserved cores
		mcbsp_util_filter_standard_pinning( machine, ret.pinning, P );
		//translate pinning if necessary
		mcbsp_util_translate_standard_pinning( machine, ret.pinning, P );
	}
	//done, return
	return ret;
}

void mcbsp_util_fatal() {
	fflush( stderr );
	exit( EXIT_FAILURE );
}

void mcbsp_util_stack_initialise( struct mcbsp_util_stack * const stack, const size_t size ) {
	stack->cap   = 16;
	stack->top   = 0;
	stack->size  = size;
	stack->array = malloc( 16 * size );
}

void mcbsp_util_stack_grow( struct mcbsp_util_stack * const stack ) {
	void * replace = malloc( 2 * stack->cap * stack->size );
	mcbsp_util_memcpy( replace, stack->array, stack->top * stack->size );
	stack->cap *= 2;
	free( stack->array );
	stack->array = replace;
}

bool mcbsp_util_stack_empty( const struct mcbsp_util_stack * const stack ) {
	return stack->top == 0;
}

void * mcbsp_util_stack_pop( struct mcbsp_util_stack * const stack ) {
	return ((char*)(stack->array)) + (stack->size * (--(stack->top)));
}

void * mcbsp_util_stack_peek( const struct mcbsp_util_stack * const stack ) {
	return ((char*)(stack->array)) + (stack->size * ( stack->top - 1 ) );
}

void mcbsp_util_stack_push( struct mcbsp_util_stack * const stack, const void * const item ) {
	if( stack->top == stack->cap )
		mcbsp_util_stack_grow( stack );

	//copy a single item into the stack
	mcbsp_util_memcpy( ((char*)(stack->array)) + stack->size * ((stack->top)++), item, stack->size );
}

void mcbsp_util_stack_destroy( struct mcbsp_util_stack * const stack ) {
	stack->cap  = 0;
	stack->top  = 0;
	stack->size = 0;
	if( stack->array != NULL ) {
		free( stack->array );
		stack->array = NULL;
	}
}

void mcbsp_util_varstack_grow( struct mcbsp_util_stack * const stack ) {
	stack->cap *= 2;
	void * replace = malloc( stack->cap );
	if( replace == NULL ) {
		fprintf( stderr, "Memory allocation failed in mcbsp_util_stack_grow!\n" );
		mcbsp_util_fatal();
	}
	mcbsp_util_memcpy( replace, stack->array, stack->top );
	free( stack->array );
	stack->array = replace;
}

void * mcbsp_util_varstack_pop( struct mcbsp_util_stack * const stack, const size_t size ) {
	stack->top -= size;
	return ((char*)(stack->array)) + stack->top;
}

void * mcbsp_util_varstack_regpop( struct mcbsp_util_stack * const stack ) {
	stack->top -= stack->size;
	return ((char*)(stack->array)) + stack->top;
}

void * mcbsp_util_varstack_peek( const struct mcbsp_util_stack * const stack, const size_t size ) {
	return ((char*)(stack->array)) + stack->top - size;
}

void * mcbsp_util_varstack_regpeek( const struct mcbsp_util_stack * const stack ) {
	return ((char*)(stack->array)) + stack->top - stack->size;
}

void mcbsp_util_varstack_push( struct mcbsp_util_stack * const stack, const void * const item, const size_t size ) {
	while( stack->top + size >= stack->cap )
		mcbsp_util_varstack_grow( stack );
	mcbsp_util_memcpy( (char*)(stack->array) + stack->top, item, size );
	stack->top += size;
}

void mcbsp_util_varstack_regpush( struct mcbsp_util_stack * const stack, const void * const item ) {
	while( stack->top + stack->size >= stack->cap )
		mcbsp_util_varstack_grow( stack );
	mcbsp_util_memcpy( (char*)(stack->array) + stack->top, item, stack->size );
	stack->top += stack->size;
}

void mcbsp_util_address_table_initialise( struct mcbsp_util_address_table * const table, const unsigned long int P ) {
	pthread_mutex_init( &(table->mutex), NULL );
	table->cap  = 16;
	table->P    =  P;
#ifdef MCBSP_ALLOW_MULTIPLE_REGS
	table->table= malloc( 16 * sizeof( struct mcbsp_util_stack * ) );
	for( unsigned long int i = 0; i < 16; ++i ) {
		table->table[ i ] = malloc( P * sizeof( struct mcbsp_util_stack ) );
		for( unsigned long int k = 0; k < table->P; ++k )
			mcbsp_util_stack_initialise( &(table->table[ i ][ k ]), sizeof(struct mcbsp_util_address_table_entry) ); 
	}
#else
	table->table= malloc( 16 * sizeof( struct mcbsp_util_address_table_entry * ) );
	for( unsigned long int i = 0; i < 16; ++i )
		table->table[ i ] = malloc( P * sizeof( struct mcbsp_util_address_table_entry ) );
#endif
}

void mcbsp_util_address_table_grow( struct mcbsp_util_address_table * const table ) {
#ifdef MCBSP_ALLOW_MULTIPLE_REGS
	struct mcbsp_util_stack ** replace = malloc( 2 * table->cap * sizeof( struct mcbsp_util_stack * ) );
	for( unsigned long int i = 0; i < 2 * table->cap; ++i ) {
		replace[ i ] = malloc( table->P * sizeof( struct mcbsp_util_stack ) );
		if( i < table->cap ) {
			for( unsigned long int k = 0; k < table->P; ++k )
				replace[ i ][ k ] = table->table[ i ][ k ];
			free( table->table[ i ] );
		} else
			for( unsigned long int k = 0; k < table->P; ++k )
				mcbsp_util_stack_initialise( &(replace[ i ][ k ]), sizeof( struct mcbsp_util_address_table_entry )  );
	}
#else
	struct mcbsp_util_address_table_entry ** replace = malloc( 2 * table->cap * sizeof( struct mcbsp_util_address_table_entry * ) );
	for( unsigned long int i = 0; i < 2 * table->cap; ++i ) {
		replace[ i ] = malloc( table->P * sizeof( struct mcbsp_util_address_table_entry ) );
		if( i < table->cap ) {
			for( unsigned long int k = 0; k < table->P; ++k )
				replace[ i ][ k ] = table->table[ i ][ k ];
			free( table->table[ i ] );
		}
	}
#endif
	free( table->table );
	table->table = replace;
	table->cap  *= 2;
}

void mcbsp_util_address_table_setsize( struct mcbsp_util_address_table * const table, const unsigned long int target_size ) {
	//if capacity is good enough, exit
	if( table->cap > target_size )
		return;

	//otherwise get lock and increase table size
	pthread_mutex_lock( &(table->mutex) );
	//check again whether another thread already
	//increased the capacity; if so, do exit
	while( table->cap <= target_size )
		mcbsp_util_address_table_grow( table );
	pthread_mutex_unlock( &(table->mutex) );
}

void mcbsp_util_address_table_destroy( struct mcbsp_util_address_table * const table ) {
	for( unsigned long int i = 0; i < table->cap; ++i ) {
#ifdef MCBSP_ALLOW_MULTIPLE_REGS
		for( unsigned long int k = 0; k < table->P; ++k )
			mcbsp_util_stack_destroy( &(table->table[ i ][ k ]) );
#endif
		free( table->table[ i ] );
	}
	free( table->table );
	pthread_mutex_destroy( &(table->mutex) );
	table->cap   = 0;
	table->P     = ULONG_MAX;
	table->table = NULL;
}

void mcbsp_util_address_table_set( struct mcbsp_util_address_table * const table, const unsigned long int key, const unsigned long int version, void * const value, const size_t size ) {
	mcbsp_util_address_table_setsize( table, key );
#ifdef MCBSP_ALLOW_MULTIPLE_REGS
	struct mcbsp_util_stack * const stack = &( table->table[ key ][ version ] );
	struct mcbsp_util_address_table_entry new_entry;
	new_entry.address = value;
	new_entry.size    = size;
	mcbsp_util_stack_push( stack, &new_entry );
#else
	struct mcbsp_util_address_table_entry * const entry = &( table->table[ key ][ version ] );
	entry->address = value;
	entry->size    = size;
#endif
}

const struct mcbsp_util_address_table_entry * mcbsp_util_address_table_get( const struct mcbsp_util_address_table * const table, const unsigned long int key, const unsigned long int version ) {
#ifdef MCBSP_ALLOW_MULTIPLE_REGS
	struct mcbsp_util_stack * const stack = &( table->table[ key ][ version ] );
	if( mcbsp_util_stack_empty( stack ) )
		return NULL;
	else {
		return (struct mcbsp_util_address_table_entry *)mcbsp_util_stack_peek( stack );
	}
#else
	return &(table->table[ key ][ version ]);
#endif
}

bool mcbsp_util_address_table_delete( struct mcbsp_util_address_table * const table, const unsigned long int key, const unsigned long int version ) {
#ifdef MCBSP_ALLOW_MULTIPLE_REGS
	struct mcbsp_util_stack * const stack = &( table->table[ key ][ version ] );
	mcbsp_util_stack_pop( stack );
	if( mcbsp_util_stack_empty( stack ) )
		return true;
	else
		return false;
#else
	struct mcbsp_util_address_table_entry * const entry = &( table->table[ key ][ version ] );
	entry->address = NULL;
	entry->size    = 0;
	return true;
#endif
}

void mcbsp_util_address_map_initialise( struct mcbsp_util_address_map * const address_map ) {
	if( address_map == NULL ) {
		fprintf( stderr, "Warning: mcbsp_util_address_map_initialise called with NULL argument.\n" );
		return;
	}
	
	address_map->cap    = 16;
	address_map->size   = 0;
	address_map->keys   = malloc( 16 * sizeof( void * ) );
	address_map->values = malloc( 16 * sizeof( size_t ) );
	
	if( address_map->keys == NULL || address_map->values == NULL ) {
		fprintf( stderr, "Error: could not allocate memory in mcbsp_util_address_map_initialise!\n" );
		mcbsp_util_fatal();
	}
}

void mcbsp_util_address_map_grow( struct mcbsp_util_address_map * const address_map ) {
	if( address_map == NULL ) {
		fprintf( stderr, "Warning: mcbsp_util_address_map_grow called with NULL argument.\n" );
		return;
	}

	//create replacement arrays
	void * * replace_k = malloc( 2 * address_map->cap * sizeof( void * ) );
	size_t * replace_v = malloc( 2 * address_map->cap * sizeof( size_t ) );

	//copy old values
	for( size_t i = 0; i < (address_map->size); ++i ) {
		replace_k[ i ] = address_map->keys[ i ];
		replace_v[ i ] = address_map->values[ i ];
	}

	//free old arrays
	free( address_map->keys   );
	free( address_map->values );

	//update struct
	address_map->keys   = replace_k;
	address_map->values = replace_v;
	address_map->cap   *= 2;
}

void mcbsp_util_address_map_destroy( struct mcbsp_util_address_map * const address_map ) {
	//delete arrays
	free( address_map->keys   );
	free( address_map->values );

	//set empty values
	address_map->cap    = 0;
	address_map->size   = 0;
	address_map->keys   = NULL;
	address_map->values = NULL;
}

size_t mcbsp_util_address_map_binsearch( const struct mcbsp_util_address_map * const address_map, const void * const key, size_t lo, size_t hi ) {

	//check boundaries
	if( address_map->keys[ lo ] == key ) return lo;
	if( address_map->keys[ hi ] == key ) return hi;

	//initialise search for insertion point
	size_t mid;
	mid  = lo + hi;
	mid /= 2;
		
	//do binary search
	do {
		//check midpoint
		if( address_map->keys[ mid ] == key )
			return mid;
		else if( key < address_map->keys[ mid ] )
			hi = mid;
		else
			lo = mid;
		//get new midpoint
		mid  = lo + hi;
		mid /= 2;
	} while( lo != mid );
	//done, return current midpoint
	return mid;
}

void mcbsp_util_address_map_insert( struct mcbsp_util_address_map * const address_map, void * const key, const size_t value ) {

	//use binary search to find the insertion point
	size_t insertionPoint;
	if( address_map->size > 0 ) {
		insertionPoint = mcbsp_util_address_map_binsearch( address_map, key, 0, address_map->size - 1 );
		if( address_map->keys[ insertionPoint ] == key ) {
			fprintf( stderr, "Warning: mcbsp_util_address_map_insert ignored as key already existed (at %lu/%lu)!\n", (unsigned long int)insertionPoint, (address_map->size-1) );
			return; //error: key was already here
		} else if( address_map->keys[ insertionPoint ] < key )
			++insertionPoint; //we need to insert at position of higher key
	} else
		insertionPoint = 0;

	//scoot over
	for( unsigned long int i = address_map->size; i > insertionPoint; --i ) {
		address_map->keys  [ i ] = address_map->keys  [ i - 1 ];
		address_map->values[ i ] = address_map->values[ i - 1 ];
	}

	//actual insert
	address_map->keys  [ insertionPoint ] = key;
	address_map->values[ insertionPoint ] = value;
	address_map->size++;

	//check capacity
	if( address_map->size == address_map->cap )
		mcbsp_util_address_map_grow( address_map );

	//done
}

void mcbsp_util_address_map_remove( struct mcbsp_util_address_map * const address_map, void * const key ) {

	if( address_map->size == 0 ) {
		fprintf( stderr, "Warning: mcbsp_util_address_map_remove ignored since map was empty!\n" );
		return; //there are no entries
	}

	//use binary search to find the entry to remove
	const size_t index = mcbsp_util_address_map_binsearch( address_map, key, 0, address_map->size - 1 );

	//check equality
	if( address_map->keys[ index ] != key ) {
		fprintf( stderr, "Warning: mcbsp_util_address_map_remove ignored since key was not found!\n" );
		return; //key not found
	}

	//scoot over to delete
	for( size_t i = index; i < address_map->size - 1; ++i ) {
		address_map->keys  [ i ] = address_map->keys  [ i + 1 ];
		address_map->values[ i ] = address_map->values[ i + 1 ];
	}

	//update size
	--(address_map->size);

	//done
}

size_t mcbsp_util_address_map_get( const struct mcbsp_util_address_map * const address_map, const void * const key ) {
	//check if we are in range
	if( address_map->size == 0 )
		return ULONG_MAX;

	//do binary search
	size_t found = mcbsp_util_address_map_binsearch( address_map, key, 0, address_map->size - 1 );

	//return if key is found
	if( address_map->keys[ found ] == key )
		return address_map->values[ found ];
	else //otherwise return max
		return SIZE_MAX;
}

void mcbsp_util_check_machine_info( struct mcbsp_util_machine_info * const machine ) {
	//sort the machine->reserved_cores list and check for duplicate entries
#ifndef NDEBUG
	machine->num_reserved_cores = mcbsp_util_sort_unique_integers( machine->reserved_cores, machine->num_reserved_cores, 0, machine->cores );
#else //performance mode
	//A quicksort usually is faster than a counting sort since typically the
	//number of reserved cores is much lower than the total number of
	//available threads.
	//Asymptotic turnover point: r log(r) > c
	const size_t turnover = (size_t)( (double)(machine->num_reserved_cores) *
		mcbsp_util_log2( machine->num_reserved_cores ) );
	if( turnover > machine->cores ) { //use counting sort
		machine->num_reserved_cores = mcbsp_util_sort_unique_integers( machine->reserved_cores, machine->num_reserved_cores, 0, machine->cores );
	} else {
		qsort( machine->reserved_cores, machine->num_reserved_cores, sizeof( size_t ), mcbsp_util_integer_compare );
	}
#endif

	//check if there are sufficient cores
	if( machine->cores <= machine->num_reserved_cores ) {
		fprintf( stderr, "Error: there are no free cores for MulticoreBSP to run on!\n" );
		mcbsp_util_fatal();
	}

	//check if we have enough threads available
	const size_t available_cores   = machine->cores - machine->num_reserved_cores;
	const size_t available_threads = available_cores * machine->threads_per_core;
	if( available_threads < machine->threads && !(machine->Tset) ) {
		fprintf( stderr, "Warning: machine information shows more threads (%lu) than #available_cores (%lu) times #threads_per_core (%lu). Attempting to continue anyway...\n",
				(unsigned long int)(machine->threads), (unsigned long int)(machine->cores), (unsigned long int)(machine->threads_per_core) );
	}

}

struct mcbsp_util_machine_partition mcbsp_util_derive_partition( const struct mcbsp_util_machine_info * const info ) {

#ifdef MCBSP_DEBUG_AFFINITY
	printf( "Given machine info: %ld threads, %ld cores, %ld threads_per_core, %ld reserved cores, %ld reserved threads_per_core.\n",
		(unsigned long int)(info->threads),
		(unsigned long int)(info->cores),
		(unsigned long int)(info->threads_per_core),
		(unsigned long int)(info->num_reserved_cores),
		(unsigned long int)(info->unused_threads_per_core)
	);
#endif

	//prepare return variable
	struct mcbsp_util_machine_partition ret;

	//start position is always 0 if we construct a partition from scratch
	ret.start = 0;

	//get maximum number of threads supported
	const size_t max_threads = info->threads;

	//machine partition is always in consecutive thread numbering. Only difference is by affinity:
	switch( info->affinity ) {
	case SCATTER:
		ret.stride = 1;
		ret.end = max_threads;
		break;
	case COMPACT:
		ret.stride = 1;
		ret.end = max_threads;
		break;
	default:
		fprintf( stderr, "Warning: unrecognised affinity in mcbsp_util_derive_partition!\n" );
	case MANUAL:
		//do not construct a partition in case of a manual or unrecognised affinity mode.
		ret.stride = 0;
		ret.end = 0;
		break;
	}

#ifdef MCBSP_DEBUG_AFFINITY
	printf( "Returning partition %ld-%ld/%ld.\n", ret.start, ret.end, ret.stride );
#endif

	//return derived partitioning
	return ret;
}

struct mcbsp_util_stack mcbsp_util_partition( const struct mcbsp_util_machine_partition subpart, const size_t P, const enum mcbsp_affinity_mode affinity ) {
	//initialise return variable
	struct mcbsp_util_stack ret;
	mcbsp_util_stack_initialise( &ret, sizeof( struct mcbsp_util_machine_partition ) );

	//floor( (end - start) / stride ) gives the number of unique available HW threads
	const size_t maxP = (subpart.end - subpart.start) / subpart.stride;

	//number of new parts
	const size_t newP = P > maxP ? maxP : P;

	//warn if we asked for too many P
	if( newP != P ) {
		fprintf( stderr, "Warning: mcbsp_util_partition is asked to partition the machine into more threads (%ld) than are available (%ld). Several BSP processes may be pinned on the same hardware thread. If this behaviour was unintended, please correct the parameters of the (nested) call(s) to bsp_begin, or supply manual affinity options.\n", P, (maxP>1?maxP:1) );
	}

	//construct each subpartition
	struct mcbsp_util_machine_partition curP;

	//loop over P
	for( size_t s = 0; s < P; ++s ) {
		//construct part s
		switch( affinity ) {
			case SCATTER:
			{
				const double fraction = ((double)maxP) / ((double)P);
				const double effectiveStride = subpart.stride * fraction;
				curP.start  = subpart.start + ((size_t)(s*effectiveStride));
				curP.stride = subpart.stride;
				curP.end    = curP.start + ((size_t)effectiveStride);
				break;
			}
			case COMPACT:
			{
				curP.start  = subpart.start  + s * subpart.stride;
				curP.stride = subpart.stride * P;
				const size_t newMaxP = (subpart.end - subpart.start) / curP.stride;
				curP.end    = curP.start + newMaxP * curP.stride;
				break;
			}
			case MANUAL:
			{
				fprintf( stderr, "Warning: mcbsp_util_partition cannot derive subpartitions when affinity is in manual mode. This function should not have been called!\n" );
				return ret;
			}
			default:
			{
				fprintf( stderr, "Error: undefined affinity in mcbsp_util_partition. Attempting to continue anyway.\n" );
				return ret;
			}
		}
		//push it on the return stack
		mcbsp_util_stack_push( &ret, &curP );
	}

	//done, return
	return ret;
}

size_t * mcbsp_util_derive_standard_pinning( const struct mcbsp_util_machine_partition * iterator, const size_t size ) {
	//allocate return array
	size_t * const ret = malloc( size * sizeof( size_t ) );
	//get a value for each entry
	for( size_t s = 0; s < size; ++s, ++iterator ) {
		//current part
		const struct mcbsp_util_machine_partition cur = *iterator;
		//get index taking the number of wraps into account
		ret[ s ] = cur.start;
	}
	//return pinning
	return ret;
}

//this has a loglinear running time
void mcbsp_util_filter_standard_pinning( const struct mcbsp_util_machine_info * const info, size_t * const array, const size_t P ) {
	//check if there are reserved cores
	if( info->num_reserved_cores == 0 )
		return; //this is not just an optimisation; the below code will fail without this check

	//filter each element in array. Note array need not be sorted.
	for( size_t i = 0; i < P; ++i ) {
		//get core number
		const size_t core = array[ i ] / info->threads_per_core; //numbering is assumed to be consecutive
		//get lower bound of current element in the reserved cores list
		const size_t offset = mcbsp_util_array_sizet_lbound( info->reserved_cores, core, 0, info->num_reserved_cores ) + 1;
		//adapt pinning accordingly
		array[ i ] += offset * info->threads_per_core;
	}
	//done
}

void mcbsp_util_translate_standard_pinning( const struct mcbsp_util_machine_info * const info, size_t * const array, const size_t P ) {
	//check if the machine already has consecutive numbering
	if( info->thread_numbering == CONSECUTIVE )
		return; //don't need to translate

	//we now expect a wrapped numbering
	if( info->thread_numbering != WRAPPED ) {
		fprintf( stderr, "Error: unrecognised thread numbering in mcbsp_util_translate_standard_pinning! Aborting...\n" );
		mcbsp_util_fatal();
	}

	//translate from a consecutive numbering to a wrapped numbering
	for( size_t i = 0; i < P; ++i ) {
		const size_t thread = array[ i ] % info->threads_per_core;
		const size_t core   = array[ i ] / info->threads_per_core;
		//translate
		array[ i ] = core + thread * info->cores;
	}
	//done
}

size_t mcbsp_util_array_sizet_lbound( const size_t * const array, const size_t key, size_t lo, size_t hi ) {

	//check boundaries
	if( array[ lo ]  > key ) return SIZE_MAX;
	if( array[ lo ] == key ) return lo;
	if( array[ hi ] == key ) return hi;

	//initialise binary search
	size_t mid = lo + hi;
	mid /= 2;

	//perform binary search
	do {
		//check midpoint
		if( array[ mid ] == key )
			return mid;
		else if( array[ mid ] > key )
			hi = mid;
		else
			lo = mid;
		//get new midpoint
		mid  = lo + hi;
		mid /= 2;
	} while( lo != mid );

	//the key did not exist; return the largest value less than the key
	return lo;
}

//does memcpy word-for-word (might be slow on architectures supporting vectorisation on larger words than chars)
static void mcbsp_util_memcpy_char( void * const restrict destination, const void * const restrict source, const size_t size ) {
	//alias input to char arrays
	const char * const restrict src = (char*)source;
	char * const restrict dest = (char*)destination;
	//do copy word-by-word
	for( size_t i = 0; i < size; ++i ) {
		dest[ i ] = src[ i ];
	}
}

//writes memcpy in terms of size_t, thus allowing more efficient vectorisation on some architectures
static void mcbsp_util_memcpy_fast( void * const restrict destination, const void * const restrict source, const size_t size ) {
	//alias input to size_t arrays
	const size_t * const restrict src = source;
	size_t * const restrict dest = destination;
	//do copy of successive data words
	size_t i = 0;
	for( ; i < size / sizeof(size_t); ++i ) {
		dest[ i ] = src[ i ];
	}
	//check if we copied all the data
	i *= sizeof(size_t);
	if( i < size ) {
		//do word-for-word copy of the remainder data
		mcbsp_util_memcpy_char( destination + i, source + i, size - i );
	}
}

//dispatches to a suitable memcpy function
void mcbsp_util_memcpy( void * const restrict destination, const void * const restrict source, const size_t size ) {
#ifdef MCBSP_USE_CUSTOM_MEMCPY
	//we do not handle memory overlap, as that is a violation of the restrict
	//contract; we will just fail hard when overlap is detected here
	assert( (destination < source && destination + size < source) ||
		(destination > source && destination + size > source) );
	//if the memory area is too small, no reason do put in much effort
	if( size < 2*sizeof(size_t) ) {
		//in this case, wide vectorisation capabilities (if available)
		//will never help accelerate the memcpy. Fall back to system
		//memcpy version
		memcpy( destination, source, size );
		//and exit
		return;
	}
	//get alignment info on destination and source pointers
	const size_t d_align = ((uintptr_t)destination) % sizeof(size_t);
	const size_t s_align = ((uintptr_t)source)      % sizeof(size_t);
	//check for matching alignments
	if( d_align == s_align ) {
		//if the matching alignment is 0
		if( d_align == 0 ) {
			//go to fast memcpy
			mcbsp_util_memcpy_fast( destination, source, size );
		} else {
			//get difference to alignment
			const size_t diff = sizeof(size_t) - d_align;
			//otherwise do word-for-word copy until aligned
			mcbsp_util_memcpy_char( destination, source, diff );
			//then proceed with fast memcpy
			mcbsp_util_memcpy_fast( destination + diff, source + diff, size - diff );
		}
	} else {
		//alignments mismatch, so just let the compiler figure it out
		mcbsp_util_memcpy_fast( destination, source, size );
	}
#else
	memcpy( destination, source, size );
#endif
}

