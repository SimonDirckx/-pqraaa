
#include "multibsp-internal.hpp"

namespace mcbsp {

	void multibsp_begin_callback() {
		MultiBSP_callable *multibsp_root_program = NULL;
		//try to get initialisation data
		const struct mcbsp_init_data * const init = static_cast< struct mcbsp_init_data * >( pthread_getspecific( mcbsp_internal_init_data ) );
		//check return value
		if( init == NULL ) {
			//check if we are in an SPMD run as a non-initialising thread
			const struct mcbsp_thread_data * const data = static_cast< struct mcbsp_thread_data * >( pthread_getspecific( mcbsp_internal_thread_data ) );
			if( data != NULL ) {
				multibsp_root_program = static_cast< MultiBSP_callable * >( data->init->bsp_program );
			} else {
				std::cerr << "Error during initialisation of the MulticoreBSP C++ wrapper; no initialisation data found (not a new BSP run), and no thread data found (not a spawned BSP process)." << std::endl;
				mcbsp_util_fatal();
			}
		} else {
			//we are the initialising thread
			multibsp_root_program = static_cast< MultiBSP_callable * >( init->bsp_program );
		}
		//set the MultiBSP program to root mode
		multibsp_root_program->_superstep = SIZE_MAX;
		//call the C library
		bsp_begin( multibsp_root_program->P );
		//first get a local child instance
		MultiBSP_callable *myInstance = multibsp_root_program->newChild();
		//run from that copy
		myInstance->spmd();
		//prevent premature destruction of our instance; other
		//processes may still communicate from local fields
		//(i.e., the bsp_direct_get, or other hp primitives.)
		bsp_sync();
		//and finally destroy that copy
		multibsp_root_program->destroyChild( myInstance );
		//exit the C SPMD section
		bsp_end();
	}

}

