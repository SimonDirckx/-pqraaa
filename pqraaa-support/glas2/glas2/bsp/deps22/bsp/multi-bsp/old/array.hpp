
#ifndef _H_MULTIBSP_ARRAY
#define _H_MULTIBSP_ARRAY

/**
 * A simple initialisable array.
 */
template< typename _T, size_t ... _tail >
class Array : public Initialisable< _tail ... > {

	protected:

		/** Fixed size. */
		size_t size;

		/** The raw array. */
		_T * array;

	public:

		/** Base constructor. */
		Array() : size( 0 ), array( NULL ) {}

		/**
		 * Default constructor.
		 * Does NOT allocate the raw array!
		 * (This should be done in a distributed way,
		 *  and not during static construction.)
		 *
		 * @param n Length of the vector.
		 */
		Array( const size_t n ) : size( n ), array( NULL ) {}

		/**
		 * Allocates the raw array locally.
		 * Should be called from within an SPMD section!
		 */
		void initialise() {
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
		void destroy() {
			//guard against uninitialised array
			if( array != NULL ) {
				//delete allocated chunk
				delete [] array;
			}
		}

		/**
		 * Return the (prospective) length of the local raw array.
		 *
		 * @return The arrray length.
		 */
		size_t length() const {
			return size;
		}

		/**
		 * Overloaded []-operator for raw array access.
		 *
		 * @param i The index of the requested element.
		 */
		_T & operator[]( const size_t i ) {
			return array[ i ];
		}
};

#endif

