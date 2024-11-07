
#ifndef _H_MULTIBSP_REPLICATED_VECTOR
#define _H_MULTIBSP_REPLICATED_VECTOR

#include "array.hpp"

/**
 * Implements a vector of a fixed size, that has its own
 * associated memory over all nodes of a Multi-BSP tree.
 * This is the leaf-case implementation; _tail is assumed
 * empty.
 *
 * @param _T The type of data this vector will store.
 * @param _tail The parameters of the Multi-BSP computer.
 */
template< typename _T, size_t ... _tail >
class Replicated_Vector_Root : public Array< _T, _tail ... > {

	public:

		/** Base constructor. */
		Replicated_Vector_Root() : Array< _T, _tail ... >() {}

		/**
		 * Default constructor.
		 * Does NOT allocate the raw array!
		 * (This should be done in a distributed way,
		 *  and not during static construction.)
		 *
		 * @param n Length of the vector.
		 */
		Replicated_Vector_Root( const size_t n ) : Array< _T, _tail ... >( n ) {}
};

template< size_t _parent_P, size_t _parent_M, typename _T, size_t ... _tail >
class Replicated_Vector : public Replicated_Vector_Root< _T, _tail... > {

	protected:

		/** Parent reference. */
		Replicated_Vector_Root< _T, _parent_P, _parent_M, _tail... > * _parent;

	public:

		Replicated_Vector() : _parent( NULL ) {}

		Replicated_Vector( const size_t n, Replicated_Vector_Root< _T, _parent_P, _parent_M, _tail... > * const _p ) :
			Replicated_Vector_Root< _T, _tail... >( n ), _parent( _p ) {}

		Replicated_Vector_Root< _T, _parent_P, _parent_M, _tail... > & parent() const {
			return *_parent;
		}

};

/**
 * Implements a vector of a fixed size, that has its own
 * associated memory over all nodes of a Multi-BSP tree.
 * This is the recursive implementation.
 *
 * @param _T The type of data this vector will store.
 * @param _tail The parameters of the Multi-BSP computer.
 */
template< typename _T, size_t _P, size_t _M, size_t ... _tail >
class Replicated_Vector_Root< _T, _P, _M, _tail... > : public Array< _T, _P, _M, _tail... > {

	protected:

		/** Child data structures. */
		Replicated_Vector< _P, _M, _T, _tail... > child_data[ _P ];

	public:

		Replicated_Vector_Root() : Array< _T, _P, _M, _tail... >() {}

		Replicated_Vector_Root( const size_t n ) : Array< _T, _P, _M, _tail... >( n ) {
			//construct static part of data structure.
			for( size_t s = 0; s < _P; ++s ) {
				child_data[ s ] = Replicated_Vector< _P, _M, _T, _tail... >( n, this );
			}
		}

		virtual Replicated_Vector< _P, _M, _T,  _tail... > & retrieve( const size_t id = bsp_pid() ) {
			return child_data[ id ];
		}

};

#endif

