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

#include "mcinternal.h"

#include <stdio.h>

void spmd( void );

void spmd( void ) {
	fprintf( stdout, "SUCCESS\n" );
}

int main() {
	if( mcbsp_internal_keys_allocated ) {
		fprintf( stderr, "FAILED \t mcbsp_internal_keys_allocated should be false!\n" );
		exit( EXIT_FAILURE );
	}
	mcbsp_internal_check_keys_allocated();
	if( !mcbsp_internal_keys_allocated ) {
		fprintf( stderr, "FAILED \t mcbsp_internal_keys_allocated should be true!\n" );
		exit( EXIT_FAILURE );
	}
	
	struct mcbsp_init_data *initialisationData = malloc( sizeof( struct mcbsp_init_data ) );
	if( initialisationData == NULL ) {
		fprintf( stderr, "FAILED \t could not allocate MulticoreBSP initialisation struct!\n" );
		exit( EXIT_FAILURE );
	}
	initialisationData->spmd = spmd;
	initialisationData->P    = 1;
	initialisationData->threadData = malloc( sizeof( void* ) );
	if( pthread_setspecific( mcbsp_internal_init_data, initialisationData ) != 0 ) {
		fprintf( stderr, "FAILED \t could not set program-specific key (pthread_setspecific)!\n" );
		exit( EXIT_FAILURE );
	}
	initialisationData = NULL;
	initialisationData = pthread_getspecific( mcbsp_internal_init_data );
	if( initialisationData == NULL ) {
		fprintf( stderr, "FAILED \t could not retrieve initialisation data (pthread_getspecific)!\n");
		exit( EXIT_FAILURE );
	}
	
	struct mcbsp_thread_data *thread_data = malloc( sizeof( struct mcbsp_thread_data ) );
	if( thread_data == NULL ) {
		fprintf( stderr, "Could not allocate local thread data!\n" );
		exit( EXIT_FAILURE );
	}
	thread_data->init   = initialisationData;
	thread_data->bsp_id = 0;
	
	mcbsp_internal_spmd( thread_data );

	mcbsp_internal_destroy_thread_data( initialisationData->threadData[ 0 ] );
	free( initialisationData->threadData );
	free( initialisationData );

	return EXIT_SUCCESS;
}

