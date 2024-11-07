
#ifndef _H_MCBSP_FFT
#define _H_MCBSP_FFT

#include <cmath>
#include <numa.h>
#include <fftw3.h>
#include <cstdlib>
#include <assert.h>
#include <iostream>
#include <stdbool.h>

#include "multibsp.hpp"
#include "mcbsp-templates.hpp"
#include "collectives/broadcast.hpp"

namespace mcbsp {
namespace numeric {

/** Common functionality for all MultiBSP_FFT_data classes. */
class FFT : public virtual mcbsp::MultiBSP_template, protected virtual mcbsp::collectives::Broadcast< fftw_plan > {

	private:

		/** Broadcast type for the FFTW plans */
		typedef mcbsp::collectives::Broadcast< fftw_plan > BCAST;

		/** Whether our current process ID is bitreversed. */
		bool bitreversed;

		/** Local vector length. */
		size_t n;

		/** Splitting point. */
		size_t k1;

		/** Leaf-level ID of this node. */
		size_t leaf_id;

		/** Mode of the FFT. */
		signed char sign;

		/** Twiddling weights. */
		double * __restrict__ tw;

		/** Buffer for use with communication during redistribution. */
		double * __restrict__ buffer;

		/** Bit-reversal permutation array of size n. */
		size_t * __restrict__ bitrev_n;

		/** Bit-reversal permutation array of size k1. */
		size_t * __restrict__ bitrev_k1;

		/** Stage 1 FFTW plan. */
		fftw_plan fftw_stage1;

		/** Stage-2 FFTW plan. */
		fftw_plan fftw_stage2;

		/** Whether our process ID is in bit-reverse. */
		bool bitrev;

		/**
		 * Initialises a bit-reversal permutation array of a given length.
		 *
		 * @param n Length of the permutation array.
		 * @param s The permutation array.
		 */
		static void bitrev_init( const size_t n, size_t * __restrict__ const s ) {
			//loop indices
			size_t j, bit;

			//constant alias: 1
			const size_t one = 1;

			//special / degenerate case
			if( n == 1 ) {
				s[ 0 ] = 0;
				return;
			}

			//general case; loop over all elements
			for( j = 0; j < n; ++j ) {
				//we should invert j and store it at s[ j ].
				//first copy j
				size_t j_copy = j;
				//init s[ j ] to 0
				s[ j ] = 0;
				//loop over all bits
				for( bit = 1; bit < n; bit <<= 1 ) {
					//shift s[ j ] to the left (to push back new bit)
					s[ j ] <<= 1;
					//copy least significant bit of j_copy into s[ j ]
					s[ j ] |= (j_copy & one);
					//shift j_copy to the right (for next bit)
					j_copy >>= 1;
				}
				//sanity check
				assert( s[ j ] < n );
			}
			//done
		}

		/** Derives the correct value for k1. */
		/** Initialises the FFTW plans. */
		void initFFTW() {
			//create a temporary array for auto-tuning
			double const * __restrict__ x = new double[ 2 * n ];
			//get the stage-1 FFT length
			const int fftw_n1 = static_cast< int >( k1 );
			//get the stage-2 FFT length
			const int fftw_n2 = static_cast< int >( n );
			//if the stage-1 is nonempty
			if( k1 > 0 ) {
				//plan this many FFTs using FFTW
				const int fftw_r1 = fftw_n2 / fftw_n1;
				fftw_stage1 = fftw_plan_many_dft( 1, &fftw_n1, fftw_r1,
						(fftw_complex*)x, &fftw_n1, 1, fftw_n1,
						(fftw_complex*)x, &fftw_n1, 1, fftw_n1,
						(sign > 0 ? -1 : 1), FFTW_PATIENT );
				//DBG
				std::cout << "Stage 1 plan created @ " << fftw_stage1 << "\n";
			}
			//if the stage-2 is nonempty
			if( k1 < n ) {
				//plan one FFT of the stage-2 length using FFTW
				fftw_stage2 = fftw_plan_dft_1d( fftw_n2,
							(fftw_complex*)x,
							(fftw_complex*)x,
							(sign > 0 ? -1 : 1), FFTW_PATIENT );
				//DBG
				std::cout << "Stage 2 plan created @ " << fftw_stage2 << "\n";
			}
			//delete our temporary array
			delete [] x;
		}

		/** Initialises a local bitrev_p array. */
		void bitrev_p_init( size_t * __restrict__ &bitrev_p ) {
			//TODO local p or global bsp_nleafs()?
			const size_t p = bsp_nprocs();
			//special case, if k1 == p
			if( k1 == p ) {
				//then bitrev_p shadows bitrev_k1
				bitrev_p = bitrev_k1;
			} else {
				//create new array
				bitrev_p = new size_t[ p ];
				//initialise values for bitreversal
				bitrev_init( p, bitrev_p );
			}
		}

		/** Initialises the weights used for twiddling in-between stage 1 and 2. */
		void twiddle_init() {
			//get the alpha corresponding to the GFFT we want to do
			const double alpha = leaf_id / static_cast< double >( k1 );
			//derive the constant part of the exponent power
			const double exp = -2.0 * M_PI * alpha / static_cast< double >( n );
			//loop over all elements of the twiddle array
			for( size_t j = 0; j < n; ++j ) {
				//set real part of the exponent
				tw[   2*j   ] = cos( j * exp );
				//set the imaginary part
				tw[ 2*j + 1 ] = sin( j * exp );
			}
		}

		/** Derives the correct value for k1. */
		void k1_init() {
			//get *total* number of processes
			const size_t  p = static_cast< size_t >(bsp_nleafs());
			//get problem size per leaf node
			const size_t np = n / p;
			//compute c such that (n/p)^c >= p
			size_t c = 1;
			for( ; c < p; c *= np ) {;}
			//store splitting stage
			k1 = n / c;
		}

		/**
		 * Code for swap-based permutations.
		 * Warning: if the permutation array is not reducible to n
		 * swaps, the permutation will fail!
		 *
		 * @param x The array which to permute.
		 * @param n The size of the array.
		 * @param s The permutation array of length n.
		 */
		static void permute( double * __restrict__ const &x, const size_t n, const size_t * __restrict__ const s ) {
			//loop index
			size_t j;
			//temporary values (imaginary, real parts)
			double tmp_i, tmp_r;
			//for each element in s
			for( j = 0; j < n; ++j ) {
				//execute the swap if this value is the lower index of the swap pair
				if( j < s[ j ] ) {
					//get old/new positions
					const size_t old_index = 2 * j;
					const size_t new_index = 2 * s[ j ];
					//cache complex value at old position
					tmp_r = x[   old_index   ];
					tmp_i = x[ old_index + 1 ];
					//put new values at old position
					x[   old_index   ] = x[   new_index   ];
					x[ old_index + 1 ] = x[ new_index + 1 ];
					//put the cached old values at the new position
					x[   new_index   ] = tmp_r;
					x[ new_index + 1 ] = tmp_i;
				}
			}
			//done
		}

		/**
		 * Redistributes the input vector x.
		 * x is assumed to be in a block distribution, and shall be
		 * redistributed to the cyclic distribution.
		 */
		void leaf_redistribute( double * __restrict__ const &x, size_t * __restrict__ const &bitrev_p ) {
			//get number of processes in this BSP run
			const unsigned int P = bsp_nprocs();
			//DBG
			/*bsp_sync();
			for( size_t k = 0; k < _P; ++k ) {
				if( k == bsp_pid() ) {
					std::cout << "pre-redistribution global array @ " << k << ": ";
					for( size_t i = 0; i < 2 * data.n; i += 2 ) {
						std::cout << *(x.pointerToElement( i )) << " ";
					}
					std::cout << std::endl;
				}
				bsp_sync();
			}*/

			//sanity check
			assert( n % P == 0 );
			//get my process ID, account for possible reversal
			const unsigned int s = bitreversed ? bitrev_p[ bsp_pid() ] : bsp_pid();
			//derive the size of each message to send to the other processes (in #doubles)
			const size_t local_size = 2 * n / P;
			//index in the buffer array
			size_t j = 0;
			//for each other process
			for( unsigned int k = 0; k < P; ++k ) {
				//remember chunk start
				double * const chunk_start = buffer + j;
				//build the message by following the
				//cyclic distribution locally
				for( size_t i = k; i < n; i += P ) {
					buffer[ j++ ] = x[   2*i   ];
					buffer[ j++ ] = x[ 2*i + 1 ]; 
				}
				//calculate offset at process k
				const size_t offset = bitreversed ? bitrev_p[ s ] * local_size : s * local_size;
				//send chunk to process k using bsp_hpput
				//for efficient buffering.
				bsp_hpput( k, chunk_start, x, offset, local_size );
			}
			//finish communication
			bsp_sync();
			//DBG
			/*bsp_sync();
			for( size_t k = 0; k < _P; ++k ) {
				if( k == bsp_pid() ) {
					std::cout << "post-redistribution global array @ " << k << ": ";
					for( size_t i = 0; i < 2 * data.n; i += 2 ) {
						std::cout << *(x.pointerToElement( i )) << " ";
					}
					std::cout << std::endl;
				}
				bsp_sync();
			}*/
		}

		/**
		 * Apply set of twiddles to the local x.
		 *
		 * @param x The input vector that needs to be twiddled.
		 */
		void twiddle( double * __restrict__ const &x ) {
			//sanity check
			assert( sign != 0 );

			//switch signs (forward if > 0, backward if < 0)
			if( sign > 0 ) {
				//loop over all local elements
				for( size_t j = 0; j < n; ++j ) {
					//get aliae of compute elements
					const double real_weight = tw[   2*j   ];
					const double imag_weight = tw[ 2*j + 1 ];
					const double real_x      =  x[   2*j   ];
					const double imag_x      =  x[ 2*j + 1 ];
					//write back twiddled factors
					x[   2*j   ] = real_weight * real_x - imag_weight * imag_x;
					x[ 2*j + 1 ] = imag_weight * real_x + real_weight * imag_x;
				}
			} else {
				//mode > 0
				//loop over all local elements
				for( size_t j = 0; j < n; ++j ) {
					//get aliae of compute elements
					const double real_weight =  tw[   2*j   ];
					//(account for sign change)
					const double imag_weight = -tw[ 2*j + 1 ];
					const double real_x      =   x[   2*j   ];
					const double imag_x      =   x[ 2*j + 1 ];
					//write back twiddled factors
					x[   2*j   ] = real_weight * real_x - imag_weight * imag_x;
					x[ 2*j + 1 ] = imag_weight * real_x + real_weight * imag_x;
				}
			}
			//done
		}

		/** Stage 1 code. */
		void stage1( double * __restrict__ const &x ) {
			//DBG
			//printVec( "Pre-stage1 local vector: " );

			//first do bit-reversal, followed by an undo last stages of the bit-reversion;
			//FFTW does not have unordered FFTs, so a call to FFTW would re-do the last
			//stages. We don't want that. First check if this entire sequence is necessary:
			if( k1 < n ) { //otherwise, k1 == n and no action is required whatsoever.
				//bit-reverse input
				permute( x, n, bitrev_n );
				//undo last stages of the bit-reversion
				for( size_t r = 0; r < n / k1 && k1 > 2; ++r ) { //if k1<2, the bit-reversal is an identity operation
					//calculate offset of this chunk
					const size_t offset = 2 * r * k1;
					permute( &(x[offset]), k1, bitrev_k1 );
				}
			} else {
				//sanity check
				assert( k1 == n );
			}

			//DBG
			std::cout << "Stage 1 executes fftw stage1 @ " << x << " with plan @ " << fftw_stage1 << "\n";

			//delegate n/k1 FFTs of length k1 to FFTW
			fftw_execute_dft( fftw_stage1, (fftw_complex*)(x), (fftw_complex*)(x) );

			//DBG
			//printVec( "Post-stage1 local vector: " );

			//done with leaf-level compute stage 1; first proceed with global redistribution.
		}

		/** Stage 2 code. */
		void stage2( double * __restrict__ const &x, size_t * __restrict__ const &bitrev_p ) {
			//do redistribution from block to cyclic
			leaf_redistribute( x, bitrev_p );

			//do a bit-reversion on x to reduce UGFFT to GFFT
			//TODO: can this be integrated into redistribute?
			permute( x, n, bitrev_n );

			//DBG
			//printVec( "Pre-twiddle local vector: " );

			//apply twiddle to reduce GFFT to FFT
			twiddle( x );

			//DBG
			//printVec( "After-twiddle local vector: " );

			//apply final stage2 FFT
			fftw_execute_dft( fftw_stage2, (fftw_complex*)(x), (fftw_complex*)(x) );

			//DBG
			//printVec( "After stage-2 local vector: " );
		}

		/** Performs top-level communication for redistribution. */
		void nested_redistribute( double * __restrict__ const &x, size_t * __restrict__ const &bitrev_p ) {
			//DBG
			/*bsp_sync();
			for( size_t k = 0; k < _P; ++k ) {
				if( k == bsp_pid() ) {
					std::cout << "pre-redistribution global array @ " << k << ": ";
					for( size_t i = 0; i < 2 * data.n; i += 2 ) {
						std::cout << x[ i ] << " ";
					}
					std::cout << std::endl;
				}
				bsp_sync();
			}*/

			//our ID
			const size_t s = bsp_pid();
			//number of leaf processes in this subtree
			const size_t subtree_P = bsp_sleafs();
			//number of processes on this level
			const size_t P = bsp_nprocs();
			//get leaf-level size
			const size_t np = n / subtree_P;
			//array of destination indices TODO move to persistent data
			size_t * __restrict__ const dest_ind = new size_t[ np ];
			//index counter
			size_t index = 0;
			//loop over blocks of size subtree_P
			for( size_t i = 0; i < n; i += subtree_P ) {
				//get lower bits
				const size_t lower = i % np;
				//get upper bits; if we are the top-level run,
				//account for bit-reversed process IDs here.
				//Then redistribution automagically fixes
				//the process IDs.
				const size_t upper_bits = i / np + s * subtree_P;
				const size_t upper = bitreversed ? bitrev_p[ upper_bits ] : upper_bits;
				//derive global index
				const size_t global = upper * np + lower;
				//derive destination process
				const size_t destproc = (global / subtree_P) % P;
				//derive destination index
				dest_ind[ index ] = ((global / subtree_P / P) * subtree_P +
							(global % subtree_P)) % n;
				//check if destination is local
				if( destproc == s && i == dest_ind[ index ] ) {
					//no action required
					continue;
				}
				//get pointer to local block
				const double * __restrict__ const block = x + 2 * i;
				//send message
				bsp_hpsend( destproc, dest_ind + index, block, 2 * subtree_P );
				//increment index
				++index;
			}
			//execute communication
			bsp_sync();
			//define pointers to tag, payloads
			const size_t * tag;
			const double * payload;
			//handle communication
			while( bsp_hpmove( (const void **)(&tag), (const void**)(&payload) ) != SIZE_MAX ) {
				//get pointer to block to be overwritten
				double * __restrict__ const block = x + 2 * *tag;
				//for each element in the payload
				for( size_t i = 0; i < 2 * subtree_P; ++i ) {
					//write element to target block
					block[ i ] = payload[ i ];
				}
			}
			//free buffer TODO move to persistent structure
			delete [] dest_ind;
			//done
			//DBG
			/*for( size_t k = 0; k < P; ++k ) {
				if( k == bsp_pid() ) {
					std::cout << "post-redistribution global array @ " << k << ": ";
					for( size_t i = 0; i < 2 * n; i += 2 ) {
						std::cout << *(x.pointerToElement( i )) << " ";
					}
					std::cout << std::endl;
				}
				bsp_sync();
			}*/
		}


	protected:

		/**
		 * Initialiser.
		 *
		 * @param _n    The local problem size.
		 * @param _sign The mode of the FFT.
		 */
		void init( const size_t _n, const signed char _sign ) {
			//initialise broadcast register
			BCAST::init();
			//initialise local fields
			if( bsp_leaf() ) {
				bitreversed = false;
				n  = _n;
				sign = _sign;
				tw = NULL;
				bitrev_n = NULL;
				bitrev_k1 = NULL;
				leaf_id = bsp_lid();
				k1_init();
				//initialise arrays
				tw        = new double[ 2 * n ];
				buffer    = new double[ 2 * n ];
				bitrev_n  = new size_t[  n ];
				bitrev_k1 = new size_t[ k1 ];
				//initialise twiddle
				twiddle_init();
				//initialise bitreversions
				bitrev_init(  n, bitrev_n  );
				bitrev_init( k1, bitrev_k1 );
				//allocate buffer
				buffer = new double[ 2 * n ];
			}
			//do FFTW autotuning if my global leaf ID is 0
			if( bsp_leaf() ) {
				if( bsp_lid() == 0 ) {
					//set libnuma to interleaved allocation
					numa_set_interleave_mask( numa_all_nodes_ptr );
					//DBG
					std::cout << bsp_lid() << " is doing initFFTW()...\n";
					initFFTW();
					//reset libnuma to local allocation
					numa_set_localalloc();
					//do broadcast; stage 1
					BCAST::leaf( &fftw_stage1, 1 );
					//stage2
					BCAST::leaf( &fftw_stage2, 1 );
					//done
				} else {
					//receive broadcasted FFTW stage 1
					fftw_stage1 = BCAST::leaf( 1 )[0];
					//stage 2
					fftw_stage2 = BCAST::leaf( 1 )[0];
				}
			} else {
				//execute broadcast code for FFTW stage 1
				BCAST::nested();
				//stage 2
				BCAST::nested();
			}
			//done
		}

		/**
		 * Destroys this leaf-level data.
		 */
		void destroy() {
			if( bsp_leaf() ) {
				//free local array, guard against empty ones
				if( tw != NULL ) {
					delete [] tw;
					tw = NULL;
				}
				if( buffer != NULL ) {
					delete [] buffer;
					buffer = NULL;
				}
				if( bitrev_n != NULL ) {
					delete [] bitrev_n;
					bitrev_n = NULL;
				}
				if( bitrev_k1 != NULL ) {
					delete [] bitrev_k1;
					bitrev_k1 = NULL;
				}
				if( buffer != NULL ) {
					delete [] buffer;
					buffer = NULL;
				}
			}
			if( bsp_leaf() && bsp_lid() == 0 ) {
				if( k1 > 0 ) {
					fftw_destroy_plan( fftw_stage1 );
				}
				if( k1 < n ) {
					fftw_destroy_plan( fftw_stage2 );
				}
			}
		}

		/** The leaf-case parallel code. */
		void leaf( double * __restrict__ const &x ) {
			//initialise bitrev_p locally
			size_t * __restrict__ bitrev_p = NULL;
			bitrev_p_init( bitrev_p );
			//do stage 1
			stage1( x );
			//perform higher-level redistributions
			bsp_up();
			//do stage 2 (includes leaf-level redistribution)
			stage2( x, bitrev_p );
			//clean up local data
			if( bitrev_p != bitrev_k1 ) {
				delete [] bitrev_p;
			}
		}

		/** The nested-case parallel code. */
		void nested( double * __restrict__ const &x ) {
			//initialise communication
			bsp_set_tagsize< size_t >();
			//initialise bitrev_p locally
			size_t * __restrict__ bitrev_p;
			bitrev_p_init( bitrev_p );
			//first do global redistributions
			bsp_up();
			//now local redistribution
			nested_redistribute( x, bitrev_p );
			//clean up local data
			if( bitrev_p != bitrev_k1 ) {
				delete [] bitrev_p;
			}
			//do lower-level redistributions plus stage 2 at leaf processes
			bsp_down();
		}
};

}; //numeric
}; //mcbsp

#endif

//TODO: template as follows:
/**
 * This class initialises data fields necessary
 * to compute an DFT of given specifications,
 * and defines functions to do this computation.
 *
 * Designed to be used within parallel FFT codes.
 *
 * This is a flat Bulk Synchronous Parallel (BSP) FFT code.
 *
 * \tparam sign The mode of the DFT. Possible values:
 *  3: Forward  mode, with scaling (1/n).
 *  2: Forward  mode, with scaling (1/sqrt(n)).
 *  1: Forward  mode, no scaling.
 * -1: Backward mode, no scaling.
 * -2: Backward mode, with scaling (1/sqrt(n)).
 * -3: Backward mode, with scaling (1/n).
 * Any other value will result in undefined
 * behaviour.
 *
 * \tparam sign      The FFT mode (forward/backward, scaling).
 * \tparam DataType  The type of floating point numbers to use (default: double).
 * \tparam LongType  The type of permutation arrays used on the input data.
 * \tparam ShortType The type of permutation arrays used on arrays of size P.
 */
//template< signed char sign, typename DataType = double, typename LongType = size_t, typename ShortType = unsigned short int >

