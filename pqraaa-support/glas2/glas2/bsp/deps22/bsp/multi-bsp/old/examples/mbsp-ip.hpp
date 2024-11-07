
#ifndef _H_MULTIBSP_IP
#define _H_MULTIBSP_IP

/**
 * Inner-product calculation according to a watered-down
 * Multi-BSP formulation.
 * (Memory sizes are not taken into account here.)
 *
 * This codes the leaf-level case.
 *
 * @tparam _tail The Multi-BSP parameters.
 */
template< size_t ... _tail, size_t _parP, size_t _parM >
class MultiBSP_IP: public mcbsp::BSP_program {

	protected:

		/** Pointer to parent class (if any). */
		const MultiBSP_IP< _P, _M,  _tail > * const parent;

		/** Problem size. */
		size_t n;

		/** Input vectors */
		Distributed_Vector< double, _tail... > &x, &y;

		virtual void spmd() {
			/** Output field. */
			double buffer[ _P ];
			//register buffer
			bsp_push_reg( buffer, _P );
			//initialise local data
			for( size_t i = 0; i < n; ++i ) {
				x[ i ] = rand() / (double)RAND_MAX;
				y[ i ] = rand() / (double)RAND_MAX;
			}
			//compute local inner product
			buffer[ bsp_pid() ] = 0.0;
			for( size_t i = 0; i < n; ++i ) {
				buffer[ bsp_pid() ] += x[ i ] * y[ i ];
			}
			//wait for sibling computations
			bsp_sync();
			//derive global inner product
			alpha = buffer[ bsp_pid() ];
			//get remote contributions and add to final result
			for( size_t s = 0; s < _P; ++s ) {
				//skip if this is my own process ID
				if( s == bsp_pid() ) continue;
				//otherwise grab data and compute inner product
				bsp_direct_get( s, buffer, s, buffer + s, 1 );
				alpha += buffer[ s ];
			}
			//done
		}

		/** Called by BSP to create a new instance for each newly created sibling SPMD process. */
		virtual BSP_program * newInstance() {
			return new MultiBSP_IP< _P, _M, _tail... >( n, parent );
		}

		/** BSP-called constructor (through newInstance). */
		MultiBSP_IP( const size_t length, 
	public:

		/** Return value. */
		double alpha;

		/** User-called constructor. */
		MultiBSP_IP(
			const size_t length,
			Distributed_Vector< double, _P, _M, _tail... > &_x,
			Distributed_Vector< double, _P, _M, _tail... > &_y ) :
			parent( NULL ), n( length ), x( _x ), y( _y ) {
			//nothing special to do here
		}

};

/** Recursive case. */
template< size_t _P, size_t _M, size_t _subP, size_t _subM, size_t ... _tail >
class MultiBSP_IP< _P, _M, _subP, _subM, _tail... >: public mcbsp::BSP_program {

	protected:

		/** Problem size. */
		size_t n;

		/** Input vectors */
		Distributed_Vector< double, _P, _M, _subP, _subM, _tail... > &x, &y;

		virtual void spmd() {
			/** Output field. */
			double buffer[ _P ];
			//register buffer
			bsp_push_reg( buffer, _P );
			//Construct nested  IP program, use nested data structures
			MultiBSP_IP< _subP, _subM, _tail... > nested( n / _P, x.retrieve(), y.retrieve() );
			//run nested program
			nested.begin( _subP );
			//get result of local nested inner product
			buffer[ bsp_pid() ] = nested.alpha;
			//wait for sibling computations
			bsp_sync();
			//derive global inner product
			alpha = buffer[ bsp_pid() ];
			//get remote contributions and add to final result
			for( size_t s = 0; s < _P; ++s ) {
				//skip if this is my own process ID
				if( s == bsp_pid() ) continue;
				//otherwise grab data and compute inner product
				bsp_direct_get( s, buffer, s, buffer + s, 1 );
				alpha += buffer[ s ];
			}
			//done
		}

		virtual BSP_program * newInstance() {
			return new MultiBSP_IP< _P, _M, _subP, _subM, _tail... >( n, x, y );
		}

	public:

		/** Return value. */
		double alpha;

		MultiBSP_IP(
			const size_t length,
			Distributed_Vector< double, _P, _M, _subP, _subM, _tail... > &_x,
			Distributed_Vector< double, _P, _M, _subP, _subM, _tail... > &_y ) :
			n( length ), x( _x ), y( _y ) {
			//nothing special to do here
		}

};


