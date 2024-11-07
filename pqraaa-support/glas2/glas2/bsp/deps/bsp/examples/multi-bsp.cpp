
#include "mcbsp.hpp"
#include "mcbsp-templates.hpp"

#include <iostream>
#include <vector>

void prettyPrintBytes( const size_t );

/**
 * A Multi-BSP datastructure that requires initialisation.
 * Recursive case.
 */
template< size_t TP, size_t TM, size_t ... Ttail >
class Initialisable {

	public:

		/** Class with virtual functions needs a virtual destructor. */
		virtual ~Initialisable() {}

		/**
		 * Only to be called from within a BSP
		 * SPMD section.
		 */
		virtual void initialise() = 0;

		virtual void destroy() = 0;

		virtual Initialisable< Ttail... > * retrievePointer() = 0;
};

/** Base case. */
template< size_t TP, size_t TM >
class Initialisable< TP, TM > {

	public:

		/** Class with virtual functions needs a virtual destructor. */
		virtual ~Initialisable() {}

		virtual void initialise() = 0;

		virtual void destroy() = 0;

		virtual void * retrievePointer() {
			std::cerr << "Error: called Initialisable::retrievePointer() on a leaf data structure!" << std::endl;
			bsp_abort( "Aborting due to illegal call to retrievePointer.\n" );
			return NULL;
		}

};

/** Recursive case. */
template< typename TT, size_t TP, size_t TM, size_t ... Ttail >
class Distributed_Vector : public Initialisable< TP, TM, Ttail... > {

	protected:

		/** Local size. */
		size_t size;

		/** Child data structures. */
		Distributed_Vector< TT, Ttail... > child_data[ TP ];

	public:

		/** Base constructor, only to be used within this class. */
		Distributed_Vector() : size( 0 ) {}

		/** Builds the static part of the data structure. */
		Distributed_Vector( const size_t n ) : size( n ) {
			for( size_t s = 0; s < TP; ++s ) {
				child_data[ s ] = Distributed_Vector< TT, Ttail... >( n / TP );
			}
		}

		/** Class with virtual functions needs a virtual destructor. */
		virtual ~Distributed_Vector() {}

		/** Reads out the data tree starting at this ndoe. */
		void readout() {
			std::cout << this << " has " << TP << " children at" << std::endl;
			for( size_t s = 0; s < TP; ++s ) {
				std::cout << "\t" << &( child_data[ s ] ) << std::endl;
			}
			for( size_t s = 0; s < TP; ++s ) {
				child_data[ s ].readout();
			}
		}

		/**
		 * Only to be called from within a BSP
		 * SPMD section.
		 */
		virtual void initialise() {
			//no action required at this level
		}

		/**
		 * Destroys the data structure at this level.
		 */
		virtual void destroy() {
			//nothing here
		}

		/** 
		 * This function can only be called from an
		 * SPMD section and on an initialised instance.
		 */
		virtual Distributed_Vector< TT, Ttail... > & retrieve() {
			return child_data[ bsp_pid() ];
		}

		/** 
		 * This function can only be called from an
		 * SPMD section and on an initialised instance.
		 */
		virtual Distributed_Vector< TT, Ttail... > * retrievePointer() {
			return &( child_data[ bsp_pid() ] );
		}

		TT * retrieveData() {
			std::cerr << "Error: trying to retrieve data from non-leaf node of this distributed data structure!" << std::endl;
		}

};

/**
 * Leaf case of the distributed vector data structure.
 * @tparam TT The data type of the vector.
 * @tparam TP The number of processes to distribute the vector over.
 * @tparam TM The available memory at this level of the Multi-BSP tree.
 */
template< typename TT, size_t TP, size_t TM >
class Distributed_Vector< TT, TP, TM > : public Initialisable< TP, TM > {

	protected:

		/** Number of elements that are stored on this level. */
		size_t size;

		/** Location of the TP arrays. */
		TT * arrays[ TP ];

	public:

		/** Base constructor. */
		Distributed_Vector() : Distributed_Vector( 0 ) {}

		/**
		 * Default constructor.
		 * @param n The total size of the array stored at this distributed vector.
		 */
		Distributed_Vector( const size_t n ) : size( n ) {
			for( size_t s = 0; s < TP; ++s ) {
				arrays[ s ] = NULL;
			}
		}

		/** Reads out the data tree starting at this ndoe. */
		void readout() {
			std::cout << this << " is a leaf node and has data fields at" << std::endl;
			for( size_t s = 0; s < TP; ++s ) {
				std::cout << "\t" << arrays[ s ] << std::endl;
			}
		}

		/**
		 * Do the actual initialisation of the data structure.
		 * This function can only be called from within an SPMD section.
		 */
		virtual void initialise() {
			//allocate size / TP (+1) elements.
			const size_t allocsize = size / TP + (size % TP > 0 ? 1 : 0);
			arrays[ bsp_pid() ] = new TT[ allocsize ];
			//check for success
			if( allocsize > 0 && arrays[ bsp_pid() ] == NULL ) {
				std::cerr << "Error: local allocation of Distributed_Vector failed!" << std::endl;
			}
		}

		/**
		 * Destroy all run-time data tied to this distributed data structure.
		 * This function can only be called from within an SPMD section.
		 */
		virtual void destroy() {
			if( arrays[ bsp_pid() ] != NULL ) {
				delete [] arrays[ bsp_pid() ];
				arrays[ bsp_pid() ] = NULL;
			}
		}

		/**
		 * @return Pointer to process-local data.
		 */
		TT * retrieveData() {
			return arrays[ bsp_pid() ];
		}

};

/**
 * This class initialises distributed data structures.
 * This class codes the non-recursive case.
 */
template< size_t TP, size_t TM, size_t ... Ttail >
class Initialiser : public mcbsp::BSP_program {

	protected:

		/** Vector of all distributed data structures to initialise. */
		std::vector< Initialisable< TP, TM > * > &toInitialise;

		/** SPMD code that does the actual initialisation. */
		virtual void spmd() {
			//iterate over all data structures
			typename std::vector< Initialisable< TP, TM > * >::iterator it = toInitialise.begin();
			for( ; it != toInitialise.end(); ++it ) {
				//call its initialisation function.
				(*it)->initialise();
			}
		}

		/** Copy constructor for use with MulticoreBSP. */
		virtual BSP_program * newInstance() {
			return new Initialiser< TP, TM >( toInitialise );
		}

	public:

		/**
		 * Defaylt constructor. 
		 * @param t The initialisable data structures that will be initialised when running this class.
		 */
		Initialiser( std::vector< Initialisable< TP, TM > * > &t ) : toInitialise( t ) {}

};

/**
 * This class initialises distributed data structures.
 * Recursive case for 4 or more template parameters.
 */
template< size_t TP, size_t TM, size_t TsubP, size_t TsubM, size_t ... Ttail >
class Initialiser< TP, TM, TsubP, TsubM, Ttail... > : public mcbsp::BSP_program {

	protected:

		/** All data structures to be initialised at this level. */
		std::vector< Initialisable< TP, TM, TsubP, TsubM, Ttail... > * > &toInitialise;

		/** SPMD code that does the actual initialisation. */
		virtual void spmd() {
			//initialise top level, build bottom-level initialisation sequence
			std::vector< Initialisable< TsubP, TsubM, Ttail... > * > bottom;
			typename std::vector< Initialisable< TP, TM, TsubP, TsubM, Ttail... > * >::iterator it = toInitialise.begin();
			for( ; it != toInitialise.end(); ++it ) {
				//call top-level init
				(*it)->initialise();
				//retrieve sub-data structure
				bottom.push_back( (*it)->retrievePointer() );
			}
			//initialise bottom level
			Initialiser< TsubP, TsubM, Ttail... > subInit = Initialiser< TsubP, TsubM, Ttail... >( bottom );
			subInit.begin( TsubP );
		}

		/** Copy constructor for use with MulticoreBSP. */
		virtual BSP_program * newInstance() {
			return new Initialiser< TP, TM, TsubP, TsubM, Ttail... >( toInitialise );
		}

	public:

		/**
		 * Defaylt constructor. 
		 * @param t The initialisable data structures that will be initialised when running this class.
		 */
		Initialiser( std::vector< Initialisable< TP, TM, TsubP, TsubM, Ttail... > * > &t ) : toInitialise( t ) {}
};

/**
 * This class destroys distributed data structures.
 * This class codes the non-recursive case.
 */
template< size_t TP, size_t TM, size_t ... Ttail >
class Destructor : public mcbsp::BSP_program {

	protected:

		/** Vector of all distributed data structures to initialise. */
		std::vector< Initialisable< TP, TM > * > &toInitialise;

		/** SPMD code that does the actual initialisation. */
		virtual void spmd() {
			//iterate over all data structures
			typename std::vector< Initialisable< TP, TM > * >::iterator it = toInitialise.begin();
			for( ; it != toInitialise.end(); ++it ) {
				//call its initialisation function.
				(*it)->destroy();
			}
		}

		/** Copy constructor for use with MulticoreBSP. */
		virtual BSP_program * newInstance() {
			return new Destructor< TP, TM >( toInitialise );
		}

	public:

		/**
		 * Defaylt constructor. 
		 * @param t The initialisable data structures that will be initialised when running this class.
		 */
		Destructor( std::vector< Initialisable< TP, TM > * > &t ) : toInitialise( t ) {}

};

/**
 * This class destroys distributed data structures.
 * Recursive case for 4 or more template parameters.
 */
template< size_t TP, size_t TM, size_t TsubP, size_t TsubM, size_t ... Ttail >
class Destructor< TP, TM, TsubP, TsubM, Ttail... > : public mcbsp::BSP_program {

	protected:

		/** All data structures to be initialised at this level. */
		std::vector< Initialisable< TP, TM, TsubP, TsubM, Ttail... > * > &toInitialise;

		/** SPMD code that does the actual initialisation. */
		virtual void spmd() {
			//initialise top level, build bottom-level initialisation sequence
			std::vector< Initialisable< TsubP, TsubM, Ttail... > * > bottom;
			typename std::vector< Initialisable< TP, TM, TsubP, TsubM, Ttail... > * >::iterator it = toInitialise.begin();
			for( ; it != toInitialise.end(); ++it ) {
				//call top-level init
				(*it)->destroy();
				//retrieve sub-data structure
				bottom.push_back( (*it)->retrievePointer() );
			}
			//initialise bottom level
			Destructor< TsubP, TsubM, Ttail... > subInit = Destructor< TsubP, TsubM, Ttail... >( bottom );
			subInit.begin( TsubP );
		}

		/** Copy constructor for use with MulticoreBSP. */
		virtual BSP_program * newInstance() {
			return new Destructor< TP, TM, TsubP, TsubM, Ttail... >( toInitialise );
		}

	public:

		/**
		 * Defaylt constructor. 
		 * @param t The initialisable data structures that will be initialised when running this class.
		 */
		Destructor( std::vector< Initialisable< TP, TM, TsubP, TsubM, Ttail... > * > &t ) : toInitialise( t ) {}
};

template< size_t TP, size_t TM, size_t ... Ttail >
class MultiBSP_IP: public mcbsp::BSP_program {

	protected:

		/** Problem size. */
		size_t n;

		/** Input vectors */
		Distributed_Vector< double, TP, TM, Ttail... > &x, &y;

		virtual void spmd() {
			//get local data
			double * const x_array = x.retrieveData();
			double * const y_array = y.retrieveData();
			//initialise local data
			for( size_t i = 0; i < n; ++i ) {
				x_array[ i ] = ((double)rand()) / ((double)RAND_MAX);
				y_array[ i ] = rand() / (double)RAND_MAX;
			}
			//compute local inner product
			alpha = 0.0;
			for( size_t i = 0; i < n; ++i ) {
				alpha += x_array[ i ] * y_array[ i ];
			}
			//register buffer
			bsp_push_reg( &alpha, 1 );
			//wait for sibling computations
			bsp_sync();
			//root process derives global inner product, through reduction
			for( unsigned int s = 1; bsp_pid() == 0 && s < TP; ++s ) {
				//local buffer
				double buffer;
				//otherwise grab data and compute inner product
 				bsp_direct_get( s, &alpha, 0, &buffer, 1 );
				alpha += buffer;
			}
			//done
		}

		virtual BSP_program * newInstance() {
			return new MultiBSP_IP< TP, TM, Ttail... >( n, x, y );
		}

	public:

		/** Return value. */
		double alpha;

		MultiBSP_IP(
			const size_t length,
			Distributed_Vector< double, TP, TM, Ttail... > &_x,
			Distributed_Vector< double, TP, TM, Ttail... > &_y ) :
			n( length ), x( _x ), y( _y ) {
			//nothing special to do here
		}

};

/** Recursive case. */
template< size_t TP, size_t TM, size_t TsubP, size_t TsubM, size_t ... Ttail >
class MultiBSP_IP< TP, TM, TsubP, TsubM, Ttail... >: public mcbsp::BSP_program {

	protected:

		/** Problem size. */
		size_t n;

		/** Input vectors */
		Distributed_Vector< double, TP, TM, TsubP, TsubM, Ttail... > &x, &y;

		virtual void spmd() {
			//Our process ID
			const unsigned int s = bsp_pid();
			//Construct nested  IP program, use nested data structures
			MultiBSP_IP< TsubP, TsubM, Ttail... > nested( n / TP, x.retrieve(), y.retrieve() );
			//register buffer
			bsp_push_reg( &(nested.alpha), 1 );
			//run nested program
			nested.begin( TsubP );
			//wait for sibling computations, make valid push_reg
			bsp_sync();
			//derive global inner product
			alpha = 0.0;
			//root process reduces remote contributions
			for( unsigned int k = 0; s == 0 && k < TP; ++k ) {
				//declare buffer
				double buffer;
				//otherwise grab data and compute inner product
				bsp_direct_get( k, &(nested.alpha), 0, &buffer, 1 );
				alpha += buffer;
			}
			//clear buffer
			bsp_pop_reg( &(nested.alpha) );
			//done
		}

		virtual BSP_program * newInstance() {
			return new MultiBSP_IP< TP, TM, TsubP, TsubM, Ttail... >( n, x, y );
		}

	public:

		/** Return value. */
		double alpha;

		MultiBSP_IP(
			const size_t length,
			Distributed_Vector< double, TP, TM, TsubP, TsubM, Ttail... > &_x,
			Distributed_Vector< double, TP, TM, TsubP, TsubM, Ttail... > &_y ) :
			n( length ), x( _x ), y( _y ) {
			//nothing special to do here
		}

};

void prettyPrintBytes( const size_t TM ) {
	if( TM >= (1l<<40) ) {
		std::cout << (TM>>40) << " TB";
	} else if( TM >= (1l<<30) ) {
		std::cout << (TM>>30) << " GB";
	} else if( TM >= (1l<<20) ) {
		std::cout << (TM>>20) << " MB";
	} else if( TM >= (1l<<10) ) {
		std::cout << (TM>>10) << " kB";
	} else {
		std::cout << TM << " B";
	}
}

template< size_t TP, size_t TM, size_t ... Ttail >
class MultiBSP_Computer_Parser {

public:

	static void print() {
		std::cout << "( " << TP << ", ";
		prettyPrintBytes( TM );
		std::cout << " ), ";
		MultiBSP_Computer_Parser< Ttail... >::print();
	}
};

template< size_t TP, size_t TM >
class MultiBSP_Computer_Parser< TP, TM > {

public:

	static void print() {
		std::cout << "( " << TP << ", ";
		prettyPrintBytes( TM );
		std::cout << " )." << std::endl;
	}
};

//Multi-BSP parameters:
//levels = 3;
//MultiBSP_P[ levels ] = { 1l, 4l, 1l }
//MultiBSP_M[ levels ] = { 1l<<15, 1l<<25, 1l<<33 }

#define MULTIBSP_TOP_P 1
#define MULTIBSP_COMPUTER MULTIBSP_TOP_P, 1l<<33, 4, 1l<<25, 1, 1l<<15

//alternative computer model
//#define MULTIBSP_TOP_P 2
//#define MULTIBSP_COMPUTER MULTIBSP_TOP_P, 1l<<33, 2, 1l<<25, 1, 1l<<15

int main() {

	std::cout << "Compiled with MultiBSP computer model:\n";
	MultiBSP_Computer_Parser< MULTIBSP_COMPUTER >::print();

	const size_t N = 1l<<18;
	Distributed_Vector< double, MULTIBSP_COMPUTER > x( N );
	Distributed_Vector< double, MULTIBSP_COMPUTER > y( N );
	
	std::vector< Initialisable< MULTIBSP_COMPUTER > * > toInit;
	toInit.push_back( &x );
	toInit.push_back( &y );

//	std::cout << "***  Pre-init read-out of the distributed vector x ***" << std::endl;
//	x.readout();

	Initialiser< MULTIBSP_COMPUTER > initialiser = Initialiser< MULTIBSP_COMPUTER >( toInit );
	initialiser.begin( MULTIBSP_TOP_P );
	
//	std::cout << std::endl << "*** Post-init read-out of the distributed vector x ***" << std::endl;
//	x.readout();

	MultiBSP_IP< MULTIBSP_COMPUTER > ip = MultiBSP_IP< MULTIBSP_COMPUTER >( N, x, y );
	ip.begin( MULTIBSP_TOP_P );

	Destructor< MULTIBSP_COMPUTER > destructor = Destructor< MULTIBSP_COMPUTER >( toInit );
	destructor.begin( MULTIBSP_TOP_P );

//	std::cout << std::endl << "*** Post-destruction read-out of the distributed vector x ***" << std::endl;
//	x.readout();

	std::cout << std::endl << "Result of inner product calculation: " << ip.alpha << std::endl;
	std::cout <<              "              Expectation value was: " << (N/4.0) << std::endl;
	return 0;
}

