
#ifndef _H_MULTIBSP_LEAF_DATA
#define _H_MULTIBSP_LEAF_DATA

/**
 * An initialisable data structure that provides each leaf in a Multi-BSP run
 * its own copy of a provided data structure type.
 *
 * @tparam The provided data structure type.
 */
template< typename _T, size_t ... _tail >
class Leaf_Data_Root : public Distributed_Vector_Root< _T, _tail... > {

	protected:

		const _T const * _template;

	public:

		/**
		 * Base constructor.
		 *
		 * Does *not* allocate any data of type _T! (Allocation should
		 * be done concurrently, and not during construction.)
		 *
		 * @param _t If not null, each leaf will have a copied instance
		 * 	     available locally. 
		 */
		Leaf_Data_Root( const _T const * _t = NULL ) : Distributed_Vector_Root< _T, _tail... >( 1 ), _template( _t ) {}

		/**
		 * Allocates local data.
		 * Should be called from within an SPMD section!
		 */
		virtual void initialise() {
			//check for template
			if( _t != NULL ) {
				//if there, just copy the template
				this->array[ 0 ] = *_t;
			} else {
				//otherwise allocate empty instance of type _T
				Distributed_Vector_Root< _T, _tail... >::initialise();
			}
		}

		/**
		 * Derefence operator.
		 * @return an instance of _T local to the current leaf processor.
		 */
		_T & operator*() {
			return *(this->array[ 0 ]);
		}

		/**
		 * Derefence operator, const version.
		 * @return an instance of _T local to the current leaf processor.
		 */
		const _T & operator*() const {
			return *(this->array[ 0 ]);
		}

		/**
		 * Dereference operator.
		 * @return a pointer to an instance of type _T which is local
		 *         to the current leaf processor.
		 */
		_T * operator->() {
			return this->array;
		}

		/**
		 * Dereference operator, const version.
		 *
		 * @return a pointer to an instance of type _T which is local
		 * to the current leaf processor.
		 */
		const _T * operator->() const {
			return this->array;
		}
};

/**
 * This class extends the previous Leaf_Data_Root class to a parent-aware
 * Leaf_Data class. This type is the one used within a Multi-BSP program.
 *
 * @tparam _parent_P The number of sibling processes; i.e., the number of
 *                   processes on the parent level.
 * @tparam _parent_M The amount of memory available on the parent level.
 * @tparam _T        The type of data this vector will store.
 * @tparam _tail     The parameters of the current Multi-BSP computer subtree.
 */
template< size_t _parent_P, size_t _parent_M, typename _T, size_t ... _tail >
class Leaf_Data : public Leaf_Data_Root< _T, _tail... > {

	protected:

		/** The parent reference. */
		Leaf_Data_Root< _T, _parent_P, _parent_M, _tail... > * parent;

	public:

		/** Base constructor (no parent reference). */
		Leaf_Data() : Distributed_Vector_Root< _T, _tail... >(), parent( NULL ) {}

		/** 
		 * Builds the static part of the data structure.
		 *
		 * @param _p A pointer to the parent data structure.
		 */
		Leaf_Data( Distributed_Vector_Root< _T, _parent_P, _parent_M, _tail... > * const _p ) :
			Leaf_Data_Root< _T, _tail... >(), _parent( _p ) {}

		/** @return A reference to the parent data structure. */
		Leaf_Data_Root< _T, _parent_P, _parent_M, _tail... > & parent() const {
			return *_parent;
		}
};

/**
 * Recursive case of Leaf_Data_Root.
 *
 * @tparam _T The type of data stored at the leaf nodes.
 * @tparam _P The top-level number of processors.
 * @tparam _M The top-level amount of available memory.
 * @tparam _tail The remaining Multi-BSP parameters.
 */
template< typename _T, size_t _P, size_t _M, size_t ... _tail >
class Leaf_Data_Root< _T, _P, _M, _tail... > : public Distributed_Vector_Root< _T, _P, _M, _tail... > {

	public:

		/** Base constructor */
		Leaf_Data_Root() : Distributed_Vector_Root< _T, _P, _M, _tail... >( MultiBSP_Computer< _P, _M, _tail >::processors() );

		/** @see Initialisable::retrieve( id ) */
		virtual Leaf_Data< _P, _M, _T, _tail... > & retrieve( const size_t id = bsp_pid() ) {
			return Distributed_Vector_Root< _T, _P, _M, _tail... >::retrieve( id );
		}
};

#endif

