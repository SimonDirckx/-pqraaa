
#ifndef _H_MULTIBSP
#define _H_MULTIBSP

#include <iostream>

extern "C" {
	#include <mcbsp.h>
}

#ifndef SIZE_MAX
 /**
  * Provide SIZE_MAX if this ANSI C99 standard macro 
  * does not make it into C++ space.
  */
 #define SIZE_MAX ((size_t)-1)
#endif

/**
 * Namespace in which the C++-style extensions for
 * MulticoreBSP for C reside.
 */
namespace mcbsp {

	/** 
	 * Shared part of any MultiBSP program. Ensures the MulticoreBSP for C
	 * back end can start SPMD sections using this program class.
	 */
	class MultiBSP_callable {

		friend void multibsp_begin_callback();

		private:

		protected:

			/** The number of threads this program should spawn. */
			size_t P;

			/** The current superstep at this level of the SPMD run. */
			size_t _superstep;

			/** The SPMD section dispatcher. */
			virtual void spmd() = 0;

			/**
			 * Method for retrieving new instances for child processes.
			 * Used to support MultiBSP_program::begin() and bsp_recurse().
			 *
			 * @return A new child instance of this MultiBSP program.
			 *
			 * @see MultiBSP_program::newChild()
			 */
			virtual MultiBSP_callable * newChild() = 0;

			/**
			 * Method for destroying instances created using newChild().
			 * Used to support MultiBSP_program::begin() and bsp_recurse().
			 *
			 * @param child The child instance to destroy.
			 *
			 * @see MultiBSP_program::destroyChild( child )
			 */
			virtual void destroyChild( MultiBSP_callable * const child ) = 0;

		public:

			/**
			 * Default constructor.
			 *
			 * @param _P The number of threads this program requires.
			 * @param _s The superstep this program will start at.
			 */
			MultiBSP_callable( size_t _P, size_t _s ) : P( _P ), _superstep( _s ) {}

	};

	/**
	 * The shared part of any implementation of a MultiBSP program. Users should
	 * derive from this MultiBSP_program class when creating their MultiBSP programs.
	 */
	template< size_t _P, size_t _M, size_t ... _tail >
	class MultiBSP_program : public MultiBSP_callable {

		private:

		protected:

			/** Number of processes to spawn on recursion. */
			size_t children;

			/** The SPMD section dispatcher. */
			virtual void spmd() {
				//check if we have control
				if( _superstep == SIZE_MAX ) {
					//yes; just go through all supersteps linearly
					for( size_t i = 0; i < max_supersteps(); ++i ) {
						superstep( i );
					}
				} else {
					//call user-defined function for current superstep
					superstep( _superstep );
					//increment superstep
					++_superstep;
					//check if out-of-bounds;
					if( _superstep == max_supersteps() ) {
						//yes, so reset counter
						_superstep = 0;
					}
				}
			}

 			/**
			 * The code for each superstep at this level.
			 * Implementation should be supplied by user.
			 *
			 * @param superstep The current superstep.
			 */
			virtual void superstep( size_t superstep ) = 0;

			/**
			 * Return the maximum number of supersteps for
			 * an SPMD run on this level.
			 * Must be supplied by user.
			 */
			virtual size_t max_supersteps() = 0;

			/**
			 * Method for retrieving new instances for sibling processes.
			 * Used to support MultiBSP_program::begin().
			 *
			 * See the BSP_program class for a detailed description.
			 */
			virtual MultiBSP_program * newChild() = 0;

			/**
			 * Method for destroying instances created using newInstance().
			 * Used to support MultiBSP_program::begin().
			 *
			 * See the BSP_program class for a detailed description.
			 *
			 * Default implementation:
			 * virtual void destroyInstance( MultiBSP_program * instance ) {
			 * 	delete instance;
			 * }
			 *
			 * @param instance The sibling program to destroy.
			 */
			virtual void destroyChild( MultiBSP_program * const instance ) {
				delete instance;
			}

		public:

			/** Default constructor. */
			MultiBSP_program() : MultiBSP_callable( 0, 0 ) {}

	};

	/** The callback function that wraps the C library around a MultiBSP::begin */
	void multibsp_begin_callback();

	/**
	 * The base template for a MultiBSP_program. User implementations should
	 * derive from this class. The code in this default template assumes
	 * only knowledge about the current node of the tree model of the target
	 * Multi-BSP computer; it thus corresponds to a leaf case. The recursive
	 * case is handled in a separate template specialisation.
	 *
	 * @tparam _P The number of processes at this level.
	 * @tparam _M The amount of available memory at this level.
	 * @tparam _tail Any remaining Multi-BSP parameters. For the implementation
	 *               here, _tail is assumed to be empty.
	 */
	template< template<size_t, size_t, size_t...> class MyProgram, size_t _P, size_t _M, size_t ... _tail >
	class MultiBSP : public MyProgram< _P, _M, _tail... > {

		private:

		protected:

		public:

			/**
			 * Starts the parallel execution of this class instance.
			 * The begin function corresponds to a top-level begin,
			 * but can be used in a nested fashion as well.
			 * If the goal is to move the execution path further
			 * down the tree of the model of the target Multi-BSP
			 * computer, however, then bsp_recurse is recommended
			 * instead.
			 */
			void begin() {
				//create a BSP-program specific initial data struct
				struct mcbsp_init_data *initialisationData = static_cast< struct mcbsp_init_data * >( malloc( sizeof( struct mcbsp_init_data ) ) );
				if( initialisationData == NULL ) {
					std::cerr << "Error: could not allocate MulticoreBSP initialisation struct!" << std::endl;
					bsp_abort( "Attemting to abort execution." );
				}
				//set values
				initialisationData->spmd        = multibsp_begin_callback;
				initialisationData->bsp_program = static_cast< void * >( this );
				initialisationData->argc        = 0;
				initialisationData->argv        = NULL;
				//continue initialisation
				bsp_init_internal( initialisationData );
				//call SPMD part
				multibsp_begin_callback();
			}

	};

	/**
	 * The base template for a MultiBSP_program. User implementations should
	 * derive from this class. The code in this specialised template assumes
	 * knowledge about the current node of the Multi-BSP tree, and of that
	 * of its child nodes as well; it thus corresponds to the case of an
	 * internal node. The leaf node case is handled in the non-specialised
	 * template definition.
	 *
	 * @tparam _P The number of processes at this level.
	 * @tparam _M The amount of available memory at this level.
	 * @tparam _subP The number of processes at the level of the child nodes.
	 * @tparam _subM The amount of memory available at the level of the child
	 * 		 nodes.
	 * @tparam _tail Any remaining Multi-BSP parameters.
	 */
	template< template< size_t, size_t, size_t, size_t, size_t ... > class MyProgram, size_t _P, size_t _M, size_t _subP, size_t _subM, size_t ... _tail >
	class MultiBSP< MyProgram, _P, _M, _subP, _subM, _tail... > : public MyProgram< _P, _M, _subP, _subM, _tail... > {

		private:

			/**
			 * The child programs.
			 */
			MultiBSP< MyProgram, _subP, _subM, _tail... > children[ _P ];

		protected:

		public:

			/**
			 * Starts the parallel execution of this class instance.
			 * The begin function corresponds to a top-level begin,
			 * but can be used in a nested fashion as well.
			 * If the goal is to move the execution path further
			 * down the tree of the model of the target Multi-BSP
			 * computer, however, then bsp_recurse is recommended
			 * instead.
			 */
			void begin() {
				//create a BSP-program specific initial data struct
				struct mcbsp_init_data *initialisationData = static_cast< struct mcbsp_init_data * >( malloc( sizeof( struct mcbsp_init_data ) ) );
				if( initialisationData == NULL ) {
					std::cerr << "Error: could not allocate MulticoreBSP initialisation struct!" << std::endl;
					bsp_abort( "Attemting to abort execution." );
				}
				//set values
				initialisationData->spmd        = multibsp_begin_callback;
				initialisationData->bsp_program = static_cast< void * >( this );
				initialisationData->argc        = 0;
				initialisationData->argv        = NULL;
				//continue initialisation
				bsp_init_internal( initialisationData );
				//call SPMD part
				multibsp_begin_callback();
			}
	};
};

template< size_t _P, size_t _M >
class BSP_Computer {

	protected:

		static void printBytes() {
			if( _M < (1l<<10 ) ) {
				std::cout << _M << " B";
			} else if( _M < (1l<<20) ) {
				std::cout << (_M>>10) << " kB";
			} else if( _M < (1l<<30) ) {
				std::cout << (_M>>20) << " MB";
			} else if( _M < (1l<<40) ) {
				std::cout << (_M>>30) << " GB";
			} else {
				std::cout << (_M>>40) << " TB";
			}
		}
};

template< size_t _P, size_t _M, size_t ... _tail >
class MultiBSP_Computer : public BSP_Computer< _P, _M > {

	public:

		static void print() {
			std::cout << "( " << _P << ", ";
			BSP_Computer< _P, _M >::printBytes();
			std::cout << " )" << std::endl;
		}

		static size_t processors() {
			return _P;
		}
};

template< size_t _P, size_t _M, size_t _subP, size_t _subM, size_t ... _tail >
class MultiBSP_Computer< _P, _M, _subP, _subM, _tail... > : public BSP_Computer< _P, _M > {

	public:

		static void print() {
			std::cout << "( " << _P << ", ";
			BSP_Computer< _P, _M >::printBytes();
			std::cout << " ), ";
			MultiBSP_Computer< _subP, _subM, _tail... >::print();
		}

		static size_t processors() {
			return _P * MultiBSP_Computer< _subP, _subM, _tail... >::processors();
		}
};

#endif

