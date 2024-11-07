
#ifndef _H_MULTIBSP_INITIALISER
#define _H_MULTIBSP_INITIALISER

#include <vector>

#include <mcbsp.hpp>

/**
 * This class initialises distributed data structures.
 * It takes pointers to multiple initialisable data structures,
 * and calls their initialise function in a single traversal of
 * the MultiBSP computer.
 * Note that the initialiser is meant to be run on-line so to
 * induce local allocation. The static part of distributed data
 * structures is handled by the Initialisable data structure
 * itself.
 *
 * This class codes the non-recursive case.
 *
 * @param _P The number of processes at this level.
 * @param _M The amount of memory available at this level.
 * @param _tail The remaining MultiBSP parameters.
 */
template< size_t _P, size_t _M, size_t ... _tail >
class Initialiser : public mcbsp::BSP_program {

	protected:

		/** Vector of all distributed data structures to initialise. */
		std::vector< Initialisable< _P, _M > * > &toInitialise;

		/** SPMD code that does the actual initialisation. */
		virtual void spmd() {
			//iterate over all data structures
			typename std::vector< Initialisable< _P, _M > * >::iterator it = toInitialise.begin();
			for( ; it != toInitialise.end(); ++it ) {
				//call its initialisation function.
				(*it)->retrieve().initialise();
			}
		}

		/** Copy constructor for use with MulticoreBSP. */
		virtual BSP_program * newInstance() {
			return new Initialiser< _P, _M >( toInitialise );
		}

	public:

		/** TODO */
		static void sequentialInit( std::vector< Initialisable< _P, _M > * > &toInitialise ) {
			typename std::vector< Initialisable< _P, _M > * >::iterator it;
			for( it = toInitialise.begin(); it != toInitialise.end(); ++it ) {
				for( size_t s = 0; s < _P; ++s ) {
					(*it)->retrieve( s ).serial_initialise();
				}
			}
		}

		/**
		 * Defaylt constructor. 
		 * @paam t The initialisable data structures that will be initialised when running this class.
		 */
		Initialiser( std::vector< Initialisable< _P, _M > * > &t ) : toInitialise( t ) {}

		void begin() {
			//first do sequential initialisations
			sequentialInit( toInitialise );
			//now defer to superclass function
			mcbsp::BSP_program::begin( _P );
		}

		void begin( size_t P ) {
			mcbsp::BSP_program::begin( P );
		}
};

/**
 * This class initialises distributed data structures.
 * Recursive case for 4 or more template parameters.
 *
 * @param _P The number of processes at this level.
 * @param _M The amount of memory at this level.
 * @param _subP The number of processes at the children level.
 * @param _subM The amount of memory each of the children have available.
 * @param _tail The remainder Multi-BSP computer parameters.
 */
template< size_t _P, size_t _M, size_t _subP, size_t _subM, size_t ... _tail >
class Initialiser< _P, _M, _subP, _subM, _tail... > : public mcbsp::BSP_program {

	protected:

		/** All data structures to be initialised at this level. */
		std::vector< Initialisable< _P, _M, _subP, _subM, _tail... > * > &toInitialise;

		/** SPMD code that does the actual initialisation. */
		virtual void spmd() {
			//initialise top level, build bottom-level initialisation sequence
			std::vector< Initialisable< _subP, _subM, _tail... > * > bottom;
			typename std::vector< Initialisable< _P, _M, _subP, _subM, _tail... > * >::iterator it;
			for( it = toInitialise.begin(); it != toInitialise.end(); ++it ) {
				//call top-level init
				(*it)->retrieve().initialise();
				//retrieve sub-data structure
				bottom.push_back( &((*it)->retrieve()) );
			}
			//initialise bottom level
			Initialiser< _subP, _subM, _tail... > subInit( bottom );
			subInit.begin( _subP );
		}

		/** Copy constructor for use with MulticoreBSP. */
		virtual BSP_program * newInstance() {
			return new Initialiser< _P, _M, _subP, _subM, _tail... >( toInitialise );
		}

	public:

		/** TODO */
		static void sequentialInit( std::vector< Initialisable< _P, _M, _subP, _subM, _tail... > * > &toInitialise ) {
			//construct next level of inits
			std::vector< Initialisable< _subP, _subM, _tail... > * > bottom;
			//prepare for iteration over initialables
			typename std::vector< Initialisable< _P, _M, _subP, _subM, _tail... > * >::iterator it;
			//for each initialisable
			for( it = toInitialise.begin(); it != toInitialise.end(); ++it ) {
				//for each process on this level
				for( size_t s = 0; s < _P; ++s ) {
					//call sequential init
					(*it)->retrieve( s ).serial_initialise();
					//construct next level
					bottom.push_back( &((*it)->retrieve( s )) );
				}
			}
			//do bottom inits
			Initialiser< _subP, _subM, _tail... >::sequentialInit( bottom );
		}

		/**
		 * Defaylt constructor. 
		 * @paam t The initialisable data structures that will be initialised when running this class.
		 */
		Initialiser( std::vector< Initialisable< _P, _M, _subP, _subM, _tail... > * > &t ) : toInitialise( t ) {}

		void begin() {
			//first do sequential initialisations
			sequentialInit( toInitialise );
			//now defer to superclass function
			mcbsp::BSP_program::begin( _P );
		}

		void begin( size_t P ) {
			mcbsp::BSP_program::begin( P );
		}
};

#endif

