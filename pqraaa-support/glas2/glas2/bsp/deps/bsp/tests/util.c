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


#include "mcutil.h"

#include <stdio.h>
#include <stdlib.h>

int main() {
	fprintf( stdout, "Warning: the exit codes for this test are switched. Failures yield a successful exit code, success results in a failure exit code.\n" );
	struct mcbsp_util_machine_info *info = mcbsp_util_createMachineInfo();
	if( info == NULL ) {
		fprintf( stderr, "FAILURE \t could not create machine info!\n" );
		exit( EXIT_SUCCESS );
	}

	info->affinity = SCATTER;
	struct mcbsp_util_pinning_info pinning = mcbsp_util_pinning( 0, info, NULL );
	if( pinning.pinning != NULL ) {
		fprintf( stderr, "FAILURE \t expected NULL value when asking for a 0-length pinning!\n" );
		exit( EXIT_SUCCESS );
	}
	free( pinning.pinning );
	mcbsp_util_stack_destroy( &(pinning.partition) );

	info->affinity = COMPACT;
	pinning = mcbsp_util_pinning( 0, info, NULL );
	if( pinning.pinning != NULL ) {
		fprintf( stderr, "FAILURE \t expected NULL value when asking for a 0-length pinning!\n" );
		exit( EXIT_SUCCESS );
	}
	mcbsp_util_stack_destroy( &(pinning.partition) );

	info->threads = info->cores = 10;
	info->affinity = SCATTER;
	pinning = mcbsp_util_pinning( 5, info, NULL );
	for( size_t i=0; i<5; ++i )
		if( pinning.pinning[ i ] != 2*i ) {
			fprintf( stderr, "FAILURE \t scatter pinning strategy does not yield expected values for P=5, on a 10-core machine!\npinning=( %ld", (unsigned long int)(pinning.pinning[ 0 ]) );
			for( size_t j=1; j<5; ++j )
				fprintf( stderr, ", %ld", (unsigned long int)(pinning.pinning[i]) );
			fprintf( stderr, " )\n" );
			exit( EXIT_SUCCESS );
		}
	free( pinning.pinning );
	mcbsp_util_stack_destroy( &(pinning.partition) );

	info->affinity  = COMPACT;
	pinning = mcbsp_util_pinning( 5, info, NULL );
	for( size_t i=0; i<5; ++i ) 
		if( pinning.pinning[ i ] != i ) {
			fprintf( stderr, "FAILURE \t compact pinning strategy does not yield expected values for P=5, on a 10-core machine!\n" );
			exit( EXIT_SUCCESS );
		}
	free( pinning.pinning );
	mcbsp_util_stack_destroy( &(pinning.partition) );

	if( info->manual_affinity != NULL ) free( info->manual_affinity );
	info->manual_affinity = malloc( info->threads * sizeof( int ) );
	info->manual_affinity[0] = 7;
	info->manual_affinity[1] = 0;
	info->manual_affinity[2] = 2;
	info->affinity = MANUAL;
	pinning = mcbsp_util_pinning( 3, info, NULL );
	if( pinning.pinning[ 0 ] != 7 || pinning.pinning[ 1 ] != 0 || pinning.pinning[ 2 ] != 2 ) {
		fprintf( stderr, "FAILURE \t manual pinning strategy does not work!\n" );
		exit( EXIT_SUCCESS );
	}
	free( pinning.pinning );
	mcbsp_util_stack_destroy( &(pinning.partition) );
	free( info->manual_affinity );
	info->manual_affinity = NULL;

	struct mcbsp_util_address_map map;

	mcbsp_util_address_map_initialise( &map );
	if( map.cap != 16 || map.size != 0 || map.keys == NULL || map.values == NULL ) {
		fprintf( stderr, "FAILURE \t initialising the mcbsp_util_address_map failed!\n" );
		exit( EXIT_SUCCESS );
	}

	mcbsp_util_address_map_grow( &map );
	if( map.cap != 32 || map.size != 0 || map.keys == NULL || map.values == NULL ) {
		fprintf( stderr, "FAILURE \t growing the mcbsp_util_address_map failed!\n" );
		exit( EXIT_SUCCESS );
	}

	mcbsp_util_address_map_insert( &map, (void*)1, 17 );
	if( map.cap != 32 || map.size != 1 || map.keys == NULL || map.values == NULL || map.keys[ 0 ] != ((void*)1) || map.values[ 0 ] != 17 ) {
		fprintf( stderr, "FAILURE \t insertion on a mcbsp_util_address_map failed!\n" );
		exit( EXIT_SUCCESS );
	}

	mcbsp_util_address_map_insert( &map, (void*)7, 13 );
	mcbsp_util_address_map_insert( &map, (void*)0, 11 );
	mcbsp_util_address_map_insert( &map, (void*)3, 30 );
	mcbsp_util_address_map_insert( &map, (void*)2, 42 );
	if( map.cap != 32 || map.size != 5 || map.keys == NULL || map.values == NULL ) {
		fprintf( stderr, "FAILURE \t repeated insertions on a mcbsp_util_address_map failed!\n" );
		exit( EXIT_SUCCESS );
	}
	for( size_t i = 1; i < map.size; ++i ) {
		if( map.keys[ i ] <= map.keys[ i - 1 ] ) {
			fprintf( stderr, "FAILURE \t insertion sort on mcbsp_util_address_map failed!\n" );
			exit( EXIT_SUCCESS );
		}
	}
	if( map.keys[ 0 ] != (void*)0 ||
		map.keys[ 1 ] != (void*)1 ||
		map.keys[ 2 ] != (void*)2 ||
		map.keys[ 3 ] != (void*)3 ||
		map.keys[ 4 ] != (void*)7 ) {
		fprintf( stderr, "FAILURE \t contents of map.keys are not as expected!\n" );
		exit( EXIT_SUCCESS );
	}
	if( map.values[ 0 ] != 11 ||
		map.values[ 1 ] != 17 ||
		map.values[ 2 ] != 42 ||
		map.values[ 3 ] != 30 ||
		map.values[ 4 ] != 13 ) {
		fprintf( stderr, "FAILURE \t contents of map.values are not as expected!\n" );
		exit( EXIT_SUCCESS );
	}
	if( mcbsp_util_address_map_get( &map, (void*)0 ) != 11 ||
		mcbsp_util_address_map_get( &map, (void*)1 ) != 17 ||
		mcbsp_util_address_map_get( &map, (void*)2 ) != 42 ||
		mcbsp_util_address_map_get( &map, (void*)3 ) != 30 ||
		mcbsp_util_address_map_get( &map, (void*)7 ) != 13 ) {
		fprintf( stderr, "FAILURE \t mcbsp_util_address_map_get does not function properly!\n" );
		exit( EXIT_SUCCESS );
	}

	mcbsp_util_address_map_remove( &map, (void*)2 );
	if( mcbsp_util_address_map_get( &map, (void*)0 ) != 11 ||
		mcbsp_util_address_map_get( &map, (void*)1 ) != 17 ||
		mcbsp_util_address_map_get( &map, (void*)2 ) != ULONG_MAX ||
		mcbsp_util_address_map_get( &map, (void*)3 ) != 30 ||
		mcbsp_util_address_map_get( &map, (void*)7 ) != 13 ) {
		fprintf( stderr, "FAILURE \t mcbsp_util_address_map_remove does not function properly!\n" );
		exit( EXIT_SUCCESS );
	}
	
	mcbsp_util_address_map_destroy( &map );
	if( map.cap != 0 || map.size != 0 || map.keys != NULL || map.values != NULL ) {
		fprintf( stderr, "FAILURE \t destroying the mcbsp_util_address_map failed!\n" );
		exit( EXIT_SUCCESS );
	}

	struct mcbsp_util_address_table table;

	mcbsp_util_address_table_initialise( &table, 5 );
	if( table.cap != 16 || table.table == NULL || table.P != 5 ) {
		fprintf( stderr, "FAILURE \t initialising the mcbsp_util_address_table failed!\n" );
		exit( EXIT_SUCCESS );
	}

	mcbsp_util_address_table_grow( &table );
	if( table.cap != 32 || table.table == NULL || table.P != 5 ) {
		fprintf( stderr, "FAILURE \t growing the mcbsp_util_address_table failed!\n" );
		exit( EXIT_SUCCESS );
	}

	mcbsp_util_address_table_set( &table, 31, 1, (void*)3,  sizeof( unsigned long int ) );
	mcbsp_util_address_table_set( &table, 64, 2, (void*)7,  sizeof( unsigned long int ) );
	mcbsp_util_address_table_set( &table, 67, 4, (void*)42, sizeof( unsigned long int ) );
#ifdef MCBSP_ALLOW_MULTIPLE_REGS
	if( ( table.cap != 128 ) ) {
		fprintf( stderr, "FAILURE \t mcbsp_util_address_table did not grow as expected!\n" );
		exit( EXIT_SUCCESS );
	}
#else
	if( ( table.cap != 128 ||
	     table.table[31][1].address != (void*)3 ||
		table.table[64][2].address != (void*)7 ||
		table.table[67][4].address != (void*)42 ) &&
		table.table[31][1].size == table.table[64][2].size &&
		table.table[64][2].size == table.table[67][4].size &&
		table.table[67][4].size == sizeof( unsigned long int ) ) {
		fprintf( stderr, "FAILURE \t mcbsp_util_address_table contents are not as expected!\n" );
		exit( EXIT_SUCCESS );
	}
#endif

	if( mcbsp_util_address_table_get( &table, 31, 1 )->address != (void*)3 ||
		mcbsp_util_address_table_get( &table, 64, 2 )->address != (void*)7 ||
		mcbsp_util_address_table_get( &table, 67, 4 )->address != (void*)42 ) {
		fprintf( stderr, "FAILURE \t mcbsp_util_address_table_get does not function properly!\n" );
		exit( EXIT_SUCCESS );
	}
	
	mcbsp_util_address_table_destroy( &table );
	if( table.cap != 0 || table.table != NULL ) {
		fprintf( stderr, "FAILURE \t destroying the mcbsp_util_address_table failed!\n" );
		exit( EXIT_SUCCESS );
	}

	struct mcbsp_util_stack stack;
	mcbsp_util_stack_initialise( &stack, sizeof( unsigned long int ) );
	if( stack.cap != 16 || stack.top != 0 || stack.size != sizeof( unsigned long int ) || stack.array == NULL ) {
		fprintf( stderr, "FAILURE \t mcbsp_util_stack did not initialise properly!\n" );
		exit( EXIT_SUCCESS );
	}

	unsigned long int toPush = 6ul;
	mcbsp_util_stack_push( &stack, &toPush );
	if( stack.cap != 16 || stack.top != 1 || stack.size != sizeof( unsigned long int ) || stack.array == NULL ||
		*((unsigned long int*)(stack.array)) != 6 ) {
		fprintf( stderr, "FAILURE \t mcbsp_util_stack_push does not function correctly!\n" );
		exit( EXIT_SUCCESS );
	}

	mcbsp_util_stack_grow( &stack );
	if( stack.cap != 32 || stack.top != 1 || stack.size != sizeof( unsigned long int ) || stack.array == NULL ||
		*((unsigned long int*)(stack.array)) != 6 ) {
		fprintf( stderr, "FAILURE \t mcbsp_util_stack_grow does not function correctly!\n" );
		exit( EXIT_SUCCESS );
	}
	
	toPush = 3ul;
	mcbsp_util_stack_push( &stack, &toPush );
	toPush = 7ul;
	for( char i = 0; i < 31; ++i )
		mcbsp_util_stack_push( &stack, &toPush );
	if( stack.cap != 64 || stack.top != 33 || stack.size != sizeof( unsigned long int ) || stack.array == NULL ||
		*((unsigned long int*)(stack.array)) != 6ul ||
		*((unsigned long int*)(((char*)stack.array)+sizeof( unsigned long int))) != 3ul ) {
		fprintf( stderr, "FAILURE \t repeated pushes on the stack fails!\n" );
		exit( EXIT_SUCCESS );
	}
	for( size_t i = 2; i < 33; ++i )
		if( *((unsigned long int*)(((char*)stack.array) + i * sizeof( unsigned long int ) )) != 7 ) {
			fprintf( stderr, "FAILURE \t repeated pushes on the stack failed!\n" );
			exit( EXIT_SUCCESS );
		}

	for( char i = 0; i < 31; ++i ) {
		if( *((unsigned long int*)(mcbsp_util_stack_peek( &stack ))) != 7ul ) {
			fprintf( stderr, "FAILURE \t mcbsp_util_stack_peek does not function properly!\n" );
			exit( EXIT_SUCCESS );
		}
		if( *((unsigned long int*)(mcbsp_util_stack_pop( &stack ))) != 7ul ) {
			fprintf( stderr, "FAILURE \t mcbsp_util_stack_pop does not function properly!\n" );
			exit( EXIT_SUCCESS );
		}
	}
	if( *((unsigned long int*)(mcbsp_util_stack_pop( &stack ))) != 3ul ||
		*((unsigned long int*)(mcbsp_util_stack_pop( &stack ))) != 6ul ||
		stack.top != 0 ) {
		fprintf( stderr, "FAILURE \t mcbsp_util_stack_pop does not return the correct values!\n" );
		exit( EXIT_SUCCESS );
	}

	mcbsp_util_stack_destroy( &stack );
	if( stack.cap != 0 || stack.top != 0 || stack.size != 0 || stack.array != NULL ) {
		fprintf( stderr, "FAILURE \t mcbsp_util_stack_destroy does not function properly!\n" );
		exit( EXIT_SUCCESS );
	}

	if( mcbsp_util_log2( 0 ) != 0 ) {
		fprintf( stderr, "FAILURE \t mcbsp_util_log2( 0 ) != 0!\n" );
		exit( EXIT_SUCCESS );
	}

	if( mcbsp_util_log2( 1 ) != 0 ) {
		fprintf( stderr, "FAILURE \t mcbsp_util_log2( 1 ) != 0!\n" );
		exit( EXIT_SUCCESS );
	}

	if( mcbsp_util_log2( 8 ) != 3 ) {
		fprintf( stderr, "FAILURE \t mcbsp_util_log2( 0 ) != 0!\n" );
		exit( EXIT_SUCCESS );
	}

	size_t unique_ints[10]={1,3,6,9,2,4,5,8,7};
	if( mcbsp_util_sort_unique_integers( unique_ints, 10, 0, 10 ) != 10 ) {
		fprintf( stderr, "FAILURE \t mcbsp_util_sort_unique_integers did not return expected number of sorted elements!\n" );
		exit( EXIT_SUCCESS );
	}
	bool sorted = true;
	for( size_t i=1; i<10; ++i ) {
		if( unique_ints[ i-1 ] >= unique_ints[ i ] ) {
			sorted = false;
			break;
		}
	}
	if( !sorted ) {
		fprintf( stderr, "FAILURE \t mcbsp_util_sort_unique_integers did not return a sorted array!\n" );
		exit( EXIT_SUCCESS );
	}

	//free up last bits of memory
	mcbsp_util_destroyMachineInfo( info );

	//Do not print "SUCCESS", the caller should verify the exit code is a failure to test this final function:
	mcbsp_util_fatal();

	return EXIT_SUCCESS;
}

