
#ifndef _H_MULTIBSP_DISTRIBUTED_VECTOR
#define _H_MULTIBSP_DISTRIBUTED_VECTOR

#include <iostream>

#include <mcbsp.hpp>

#include "initialisable.hpp"

/**
 * Implements a vector distributed over a the hierarchy
 * of a Multi-BSP computer. This class represents the
 * base type on a single compute node; that is, a raw
 * array of the requested type.
 *
 * @tparam _T The type of data this vector will store.
 * @tparam _tail The parameters of the Multi-BSP computer.
 */
template< typename _T, size_t ... _tail >
class Distributed_Vector_Root : public Initialisable< _tail ... > {

	protected:

		/** Local size. */
		size_t size;

		/** The raw array. */
		_T * array;

	public:

		/** Base constructor. */
		Distributed_Vector_Root() : size( 0 ), array( NULL ) {}

		/** 
		 * Default constructor.
		 * Does NOT allocate the raw array!
		 * (This should be done in a distributed way,
		 *  and not during static construction.)
		 *
		 * @param n Length of the vector.
		 */
		Distributed_Vector_Root( const size_t n ) : size( n ), array( NULL ) {}

		/**
		 * Allocates the raw array locally.
		 * Should be called from within an SPMD section!
		 */
		virtual void initialise() {
			//guard against empty alloc
			if( size > 0 ) {
				//do alloc
				array = new _T[ size ];
				//check alloc
				if( array == NULL ) {
					bsp_abort( "Error: could not allocate raw array of requested length %ld!\n", size );
				}
			}
		}

		/**
		 * Frees the locally allocated raw array.
		 * Should be called from within an SPMD section!
		 */
		virtual void destroy() {
			//guard against uninitialised array
			if( array != NULL ) {
				//delete allocated chunk
				delete [] array;
			}
		}

		/**
		 * Return the (prospective) length of the
		 * local raw array.
		 * 
		 * @return The arrray length.
		 */
		size_t length() const {
			return size;
		}

		/** Reads out the data tree starting at this node. */
		void readout() const {
			std::cout << this << " wraps a raw array of length " << size << " at address " << array << std::endl;
		}

		/**
		 * Overloaded []-operator for raw array access.
		 *
		 * @param i The index of the requested element.
		 */
		_T & operator[]( const size_t i ) {
			return array[ i ];
		}

		/**
		 * Gets a pointer to an element that is distributed
		 * somewhere in this distributed (sub)vector.
		 *
		 * @param vector The given distributed (sub)vector.
		 * @param i      The global index in the distributed (sub)vector.
		 */
		_T * pointerToElement( const size_t i ) {
			//sanity check
			assert( i < size );
			//return requested pointer
			return array + i;
		}
};

/**
 * Implements a vector distributed over a the hierarchy
 * of a Multi-BSP computer. This class extends the
 * previous Distributed_Vector_Root class to a parent-
 * aware Distributed_Vector.
 * This type is the one used within a Multi-BSP program.
 *
 * @tparam _parent_P The number of sibling processes; i.e., the number of processes on the parent level.
 * @tparam _parent_M The amount of memory available on the parent level.
 * @tparam _T The type of data this vector will store.
 * @tparam _tail The parameters of the current Multi-BSP computer subtree.
 */
template< size_t _parent_P, size_t _parent_M, typename _T, size_t ... _tail >
class Distributed_Vector : public Distributed_Vector_Root< _T, _tail... > {

	protected:

		/** Parent reference. */
		Distributed_Vector_Root< _T, _parent_P, _parent_M, _tail... > * _parent;

	public:

		/** No parent reference in default constructor. */
		Distributed_Vector() : _parent( NULL ) {}

		/**
		 * Builds the static part of the data structure.
		 *
		 * @param _n The size of the array at this level.
		 * @param _p Pointer to the parent data structure.
		 */
		Distributed_Vector( const size_t _n, Distributed_Vector_Root< _T, _parent_P, _parent_M, _tail... > * const _p ) :
			Distributed_Vector_Root< _T, _tail... >( _n ), _parent( _p ) {}

		/** Gets a reference to the parent. */
		Distributed_Vector_Root< _T, _parent_P, _parent_M, _tail... > & parent() const {
			return *_parent;
		}
};

/**
 * Recursive case. Implements a vector distributed over
 * a the hierarchy of a Multi-BSP computer.
 *
 * @tparam _T The type of data this vector will store.
 * @tparam _P The top-level number of processors.
 * @tparam _M The top-level amount of available memory.
 * @tparam _tail The number of processors and available
 *              memory of the remaining levels of the
 *              Multi-BSP hierarchy.
 */
template< typename _T, size_t _P, size_t _M, size_t ... _tail >
class Distributed_Vector_Root< _T, _P, _M, _tail... > : public Initialisable< _P, _M, _tail... > {

	protected:

		/** Local size. */
		size_t size;

		/** Child data structures. */
		Distributed_Vector< _P, _M, _T, _tail... > child_data[ _P ];

	public:

		/** Base constructor, only to be used within this class. */
		Distributed_Vector_Root() : size( 0 ) {}

		/** Builds the static part of the data structure. */
		Distributed_Vector_Root( const size_t n ) : size( n ) {
			//construct static part of data structure.
			for( size_t s = 0; s < _P; ++s ) {
				//account for change in size of lower-level array
				child_data[ s ] = Distributed_Vector< _P, _M, _T, _tail... >( n / _P, this );
			}
		}

		/** Reads out the data tree starting at this ndoe. */
		void readout() const {
			std::cout << this << " has " << _P << " children at" << std::endl;
			for( size_t s = 0; s < _P; ++s ) {
				std::cout << "\t" << &( child_data[ s ] ) << std::endl;
			}
			for( size_t s = 0; s < _P; ++s ) {
				child_data[ s ].readout();
			}
		}

		/**
		 * Initialises the distributed data structure
		 * at this level.
		 * Must be called from within an SPMD section.
		 */
		virtual void initialise() {
			//no action required at this level
		}

		/**
		 * Destroys the data structure at this level.
		 * Must be called from within an SPMD section.
		 */
		virtual void destroy() {
			//nothing here
		}

		/**
		 * Return the (prospective) length of the
		 * local raw array.
		 * 
		 * @return The arrray length.
		 */
		size_t length() const {
			return size;
		}

		/**
		 * @see Initialisable::retrieve( id )
		 */
		virtual Distributed_Vector< _P, _M, _T, _tail... > & retrieve( const size_t id = bsp_pid() ) {
			return child_data[ id ];
		}

		/**
		 * Gets a pointer to an element that is distributed
		 * somewhere in this distributed (sub)vector.
		 *
		 * @param i The global index in the distributed (sub)vector.
		 */
		_T * pointerToElement( const size_t i ) {
			//sanity check
			assert( i < size );
			//chunk size
			const size_t chunk = size / _P;
			//determine child
			const size_t child = i / chunk;
			//recurse there
			return retrieve( child ).pointerToElement( i % chunk );
		}
};

#endif

