
#include <vector>
#include <assert.h>

#include <multibsp.hpp>
#include <mcbsp-templates.hpp>
#include <distributed-vector.hpp>
#include <replicated-vector.hpp>
#include <initialiser.hpp>
#include <destructor.hpp>

#ifndef _H_MULTIBSP_IP
#define _H_MULTIBSP_IP

template< size_t _P >
class BSPIP_common {

	protected:

		/** Field where the local result will be stored. */
		double local_result;

		/**
		 * Communication pattern that is used for all-
		 * to-all transmission of the local results,
		 * and the summation of each of those local
		 * results into the final result value alpha.
		 *
		 * WARNING: assumes local_result is registered
		 * 	    for communication!
		 */
		void allreduce( double &alpha ) {
			//buffer for received data
			double recv_buffer;
			//derive global inner product
			alpha = local_result;
			//get remote contributions and add to final result
			for( size_t s = 0; s < _P; ++s ) {
				//skip if this is my own process ID
				if( s == bsp_pid() ) continue;
				//otherwise grab data
				bsp_direct_get( s, &local_result, 0, &recv_buffer );
				//and compute inner product
				alpha += recv_buffer;
			}
		}
};

//Note: the following macro would expand to the same two lines:
//MULTIBSP_LEAF_CLASS( MultiBSP_IP ) {
template< size_t _P, size_t _M, size_t ... _tail >
class MultiBSP_IP : public mcbsp::MultiBSP_program< _P, _M, _tail... >, public BSPIP_common< _P > {
	
		MULTIBSP_LEAF_HEADER

	protected:

		Distributed_Vector< _P, _M, double, _tail... > &x, &y;
		Replicated_Vector< _P, _M, double, _tail... > &alpha;

		virtual MultiBSP_IP< _P, _M, _tail... > * newInstance() {
			return new MultiBSP_IP< _P, _M, _tail... >( x.parent(), y.parent(), alpha.parent(), bsp_pid() );
		}

		virtual void superstep( size_t superstep ) {
			if( superstep == 0 ) {
				//get length
				const size_t n = x.length();
				//register buffer
				bsp_push_reg( &(this->local_result) );
				//initialise local data
				for( size_t i = 0; i < n; ++i ) {
					x[ i ] = rand() / (double)RAND_MAX;
					y[ i ] = rand() / (double)RAND_MAX;
				}
				//compute local inner product
				this->local_result = 0.0;
				for( size_t i = 0; i < n; ++i ) {
					 this->local_result += x[ i ] * y[ i ];
				}
				//wait for sibling computations
				bsp_sync();
				//do all-to-all and reduction
				this->allreduce( alpha[ 0 ] );
				//done
			}
		}

		//there is only one superstep here
		virtual size_t max_supersteps() { return 1; }

	public:

		MultiBSP_IP( Distributed_Vector_Root< double, _P, _M, _tail... > &_x,
				Distributed_Vector_Root< double, _P, _M, _tail... > &_y,
				Replicated_Vector_Root< double, _P, _M, _tail... > &_alpha,
				const size_t id = 0 ) :
			x( _x.retrieve( id ) ), y( _y.retrieve( id ) ), alpha( _alpha.retrieve( id ) ) {}
};


//Note: the following macro would expand to the same two lines:
//       MULTIBSP_INTERNAL_CLASS( MultiBSP_IP ) {
template< size_t _P, size_t _M, size_t _subP, size_t _subM, size_t ... _tail >
class MultiBSP_IP< _P, _M, _subP, _subM, _tail... > : public mcbsp::MultiBSP_program< _P, _M, _subP, _subM, _tail... >, public BSPIP_common< _P > {

	MULTIBSP_INTERNAL_HEADER

	protected:

		Distributed_Vector< _P, _M, double, _subP, _subM, _tail... > &x, &y;
		 Replicated_Vector< _P, _M, double, _subP, _subM, _tail... > &alpha;

		virtual void spmd() {
			//register buffer
			bsp_push_reg( &(this->local_result) );
			//recuse to compute first
			bsp_recurse();
			//retrieve solution
			this->local_result = alpha.retrieve( 0 )[ 0 ];
			//allow all children subprograms to finish
			bsp_sync();
			//do allreduce
			this->allreduce( alpha[ 0 ] );
		}

		virtual void superstep( size_t superstep ) {
			if( superstep == 0 ) {
				//register buffer
				bsp_push_reg( &(this->local_result) );
				//recuse to compute first
				bsp_recurse();
				//retrieve solution
				this->local_result = alpha.retrieve( 0 )[ 0 ];
			} else {
				assert( superstep == 1 );
				//do allreduce
				this->allreduce( alpha[ 0 ] );
			}
		}

		virtual size_t max_supersteps() { return 2; }		

		virtual MultiBSP_IP< _P, _M, _subP, _subM, _tail... > * newInstance() {
			return new MultiBSP_IP< _P, _M, _subP, _subM, _tail... >( x.parent(), y.parent(), alpha.parent(), bsp_pid() );
		}

		virtual MultiBSP_IP< _subP, _subM, _tail... > * newChild( const size_t s ) {
			return new MultiBSP_IP< _subP, _subM, _tail... >( x, y, alpha, s );
		}

	public:

		/**
		 * Base constructor. Takes the parent data structure as input, and automatically
		 * retrieves the local data structures according to a given process ID.
		 * 
		 * @param _x The left input vector.
		 * @param _y The right input vector.
		 * @param id The ID of the thread calling the constructor. (Defaults to ID 0.)
		 */
		MultiBSP_IP( Distributed_Vector_Root< double, _P, _M, _subP, _subM, _tail... > &_x,
				Distributed_Vector_Root< double, _P, _M, _subP, _subM, _tail... > &_y,
				 Replicated_Vector_Root< double, _P, _M, _subP, _subM, _tail... > &_alpha,
				const size_t id = 0 ) :
			x( _x.retrieve( id ) ), y( _y.retrieve( id ) ), alpha( _alpha.retrieve( id ) ) {}

};

#endif

