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

#include "mcbsp-internal.hpp"

namespace mcbsp {

	void BSP_program::begin( const MCBSP_PROCESSOR_INDEX_DATATYPE _P ) {
		//set number of threads to use
		P = static_cast< size_t >( _P );
		//create a BSP-program specific initial data struct
		struct mcbsp_init_data *initialisationData = static_cast< struct mcbsp_init_data * >( malloc( sizeof( struct mcbsp_init_data ) ) );
		if( initialisationData == NULL ) {
			fprintf( stderr, "Error: could not allocate MulticoreBSP initialisation struct!\n" );
			mcbsp_util_fatal();
		}
		//set values
		initialisationData->spmd 	= &mcbsp_cpp_callback;
		initialisationData->bsp_program = static_cast< void * >( this );
		initialisationData->argc 	= 0;
		initialisationData->argv 	= NULL;

		//continue initialisation
		bsp_init_internal( initialisationData );

		//call SPMD part
		mcbsp_cpp_callback();
	}

	void BSP_program::destroyInstance( BSP_program * const instance ) {
		delete instance;
	}

	BSP_program::~BSP_program() {}

	void mcbsp_cpp_callback() {
		BSP_program *bsp_program = NULL;
		//try to get initialisation data
		const struct mcbsp_init_data * const init = static_cast< struct mcbsp_init_data * >( pthread_getspecific( mcbsp_internal_init_data ) );
		//check return value
		if( init == NULL ) {
			//check if we are in an SPMD run as a non-initialising thread
			const struct mcbsp_thread_data * const data = static_cast< struct mcbsp_thread_data * >( pthread_getspecific( mcbsp_internal_thread_data ) );
			if( data != NULL ) {
				bsp_program = static_cast< BSP_program * >( data->init->bsp_program );
			} else {
				std::cerr << "Error during initialisation of the MulticoreBSP C++ wrapper; no initialisation data found (not a new BSP run), and no thread data found (not a spawned BSP process)." << std::endl;
				mcbsp_util_fatal();
			}
		} else {
			//we are the initialising thread
			bsp_program = static_cast< BSP_program * >( init->bsp_program );
		}
		//retieve P
		const size_t P = bsp_program->P;

		//call the C library
		bsp_begin( static_cast< MCBSP_PROCESSOR_INDEX_DATATYPE >( P ) );
		//if we are the initialising thread, call the original class instance
		if( bsp_pid() == 0 ) {
			bsp_program->spmd();
			//matching sync for the else-case
			bsp_sync();
		} else {
			//otherwise, first get a copy
			BSP_program *myInstance = bsp_program->newInstance();
			//run from that copy
			myInstance->spmd();
			//prevent premature destruction of our instance; other
			//processes may still communicate from local fields
			//(i.e., the bsp_direct_get, or other hp primitives.)
			bsp_sync();
			//and finally destroy that copy
			bsp_program->destroyInstance( myInstance );
		}
		//exit the C SPMD section
		bsp_end();
	}

}

