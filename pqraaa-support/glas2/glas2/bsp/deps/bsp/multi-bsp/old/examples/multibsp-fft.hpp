
#include <math.h>
#include <assert.h>

#include <fftw3.h>

#include <multibsp.hpp>
#include <mcbsp-templates.hpp>
#include <distributed-vector.hpp>
#include <initialiser.hpp>
#include <destructor.hpp>

/** Common functionality for all MultiBSP_FFT_data classes. */
class MultiBSP_FFT_data_common {

	protected:

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

	public:

		/** Default constructor. */
		MultiBSP_FFT_data_common() : bitrev_p( NULL ) {}

		/** Bit-reversal permutation array of size p. */
		size_t * __restrict__ bitrev_p;
};

/**
 * Data corresponding to the MultiBSP FFT algorithm.
 * This class codes the leaf data.
 *
 * @tparam _tail The remaining Multi-BSP parameters
 */
template< size_t ... _tail >
class MultiBSP_FFT_data : public Initialisable< _tail... >, public MultiBSP_FFT_data_common {

	protected:

		/** Derives the correct value for k1. */
		void k1_init() {
			//store splitting stage (degenerate case)
			k1 = n;
		}

		/**
		 * Direct constructor.
		 *
		 * @param _n    The local problem size.
		 * @param _k1   The splitting stage.
		 * @param _sign The mode of the FFT.
		 * @param _s    The process ID corresponding to this leaf node.
		 * @param _leaf The parent's leaf ID.
		 * @param _PIDs The PIDs of all parents.
		 * @param _Ps   The number of processes on each parent level.
		 */
		MultiBSP_FFT_data( const size_t _n, const size_t _k1, const signed char _sign,
					const size_t _s, const size_t _leaf, const std::vector< size_t > _PIDs,
					const std::vector< size_t > _Ps ) :
			n( _n ), k1( _k1 ), superstep( 0 ), sign( _sign ), PIDs( _PIDs ), Ps( _Ps ),
			tw( NULL ), bitrev_n( NULL ), bitrev_k1( NULL ) {
			//record current PID
			PIDs.push_back( _s );
			//derive leaf ID
			leaf_id = _leaf + _s;
		}

		/**
		 * Default constructor. Initialises with meaningless data.
		 * Only to be used internally as a placeholder.
		 */
		MultiBSP_FFT_data() : n( 0 ), k1( 0 ), superstep( 0 ), leaf_id( 0 ),  sign( 0 ), tw( NULL ), bitrev_n( NULL ), bitrev_k1( NULL ) {}

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
			}
			//if the stage-2 is nonempty
			if( k1 < n ) {
				//plan one FFT of the stage-2 length using FFTW
				fftw_stage2 = fftw_plan_dft_1d( fftw_n2,
							(fftw_complex*)x,
							(fftw_complex*)x,
							(sign > 0 ? -1 : 1), FFTW_PATIENT );
			}
			//delete our temporary array
			delete [] x;
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

	public:

		/** Local vector length. */
		size_t n;

		/** Splitting point. */
		size_t k1;

		/** Current superstep. */
		size_t superstep;

		/** Leaf-level ID of this node. */
		size_t leaf_id;

		/** Mode of the FFT. */
		signed char sign;

		/** Process IDs. */
		std::vector< size_t > PIDs;

		/** Number of processes on the parent level. */
		std::vector< size_t > Ps;

		/** Twiddling weights. */
		double * __restrict__ tw;

		/** Buffer for use with communication. */
		double * __restrict__ buffer;

		/** Bit-reversal permutation array of size n. */
		size_t * __restrict__ bitrev_n;

		/** Bit-reversal permutation array of size k1. */
		size_t * __restrict__ bitrev_k1;

		/** Stage 1 FFTW plan. */
		fftw_plan fftw_stage1;

		/** Stage-2 FFTW plan. */
		fftw_plan fftw_stage2;

		/**
		 * Base constructor. Automatically derives k1.
		 * Initialises with process ID assumed to be 0.
		 *
		 * @param _n    The local problem size.
		 * @param _sign The mode of the FFT.
		 * @param _x    The local input vector.
		 */
		MultiBSP_FFT_data( const size_t _n, const signed char _sign ) :
			n( _n ), superstep( 0 ), leaf_id( 0 ), sign( _sign ), tw( NULL ), buffer( NULL ), bitrev_n( NULL ), bitrev_k1( NULL ) {
			k1_init();
		}

		/**
		 * Initialises this leaf level of the FFTW data.
		 * Should be called from within an SPMD section.
		 */
		virtual void initialise() {
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

			//only init bitrev_p if we are the top-level run
			if( superstep == 0 ) {
				//special case; if k1 == p
				const size_t p = bsp_nprocs();
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
		}

		/**
		 * Does FFTW auto-tuning.
		 * Done here since FFTW auto-tuning is not thread-safe.
		 */
		virtual void serial_initialise() {
			//initialise FFTW plans
			initFFTW();
		}

		/**
		 * Destroys this leaf-level data.
		 */
		virtual void destroy() {
			//destroy FFTW data (assume this is thread-safe)
			//fftw_destroy_plan( fftw_stage1 );
			//fftw_destroy_plan( fftw_stage2 );
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
			if( bitrev_p != NULL ) {
				//check if we are not in special case
				if( bitrev_p != bitrev_k1 ) {
					//delete if we have our own array
					delete [] bitrev_p;
				}
				//in any case, set to NULL
				bitrev_p = NULL;
			}
			if( bitrev_k1 != NULL ) {
				delete [] bitrev_k1;
				bitrev_k1 = NULL;
			}
		}

		/**
		 * Destroys the FFTW auto-tuning plans.
		 * Done here since FFTW's destroy functions are not thread-safe.
		 */
		virtual void serial_destroy() {
			//do destruction here, as FFTW is also not thread-safe here
			if( k1 > 0 ) {
				fftw_destroy_plan( fftw_stage1 );
			}
			if( k1 < n ) {
				fftw_destroy_plan( fftw_stage2 );
			}
		}
};

/**
 * Data corresponding to the MultiBSP FFT algorithm.
 * This class codes the internal, non-root, data.
 *
 * @tparam _P    The number of processes at the parent level.
 * @tparam _M    The amount of memory at the parent level.
 * @tparam _tail The remaining Multi-BSP parameters
 */
template< size_t _P, size_t _M, size_t ... _tail >
class MultiBSP_FFT_internal : public MultiBSP_FFT_data< _tail... > {

	friend class MultiBSP_FFT_data< _P, _M, _tail... >;

	protected:

		/** Parent reference. */
		MultiBSP_FFT_data< _P, _M, _tail... > * parent_pointer;

		/**
		 * Builds the static part of the data structure.
		 * Direct constructor.
		 *
		 * @param _n    The local problem size.
		 * @param _k1	The splitting stage.
		 * @param _sign The mode of the FFT.
		 * @param _s    The process ID corresponding to this leaf node.
		 * @param _leaf The parent's leaf ID.
		 * @param _PIDs The PIDs of all parents.
		 * @param _Ps   The number of processes on the parent level.
		 */
		MultiBSP_FFT_internal( const size_t _n, const size_t _k1, const signed char _sign,
					const size_t _s, const size_t _leaf,
					const std::vector< size_t > _PIDs,
					const std::vector< size_t > _Ps,
					MultiBSP_FFT_data< _P, _M, _tail... > *_p ) :
			MultiBSP_FFT_data< _tail... >( _n, _k1, _sign, _s, _leaf, _PIDs, _Ps ), parent_pointer( _p ) {}

	public:

		/** Default constructor. */
		MultiBSP_FFT_internal() : parent_pointer( NULL ) {}

		/**
		 * Builds the static part of the data structure.
		 * Automatically derives k1.
		 *
		 * @param _n    The local problem size.
		 * @param _sign The mode of the FFT.
		 * @param _x    The local input vector.
		 */
		MultiBSP_FFT_internal( const size_t _n, const signed char _sign,
					MultiBSP_FFT_data< _P, _M, _tail... > *_p ) :
			MultiBSP_FFT_data< _tail... >( _n, _sign ), parent_pointer( _p ) {}

		/** Gets a reference to the parent. */
		MultiBSP_FFT_data< _P, _M, _tail... > & parent() const {
			return *parent_pointer;
		}
};

/**
 * Data corresponding to the MultiBSP FFT algorithm.
 * This class codes the internal data.
 *
 * @tparam _P    The number of child processes.
 * @tparam _M    The amount of memory available at each child.
 * @tparam _tail The remaining Multi-BSP parameters
 */
template< size_t _P, size_t _M, size_t ... _tail >
class MultiBSP_FFT_data< _P, _M, _tail... > : public Initialisable< _P, _M, _tail... >, public MultiBSP_FFT_data_common {

	protected:

		/** Child data structures. */
		MultiBSP_FFT_internal< _P, _M, _tail... > children[ _P ];

		/** Derives the correct value for k1. */
		void k1_init() {
			//get *total* number of processes
			const size_t  p = MultiBSP_Computer< _P, _M, _tail... >::processors();
			//get problem size per leaf node
			const size_t np = n / p;
			//compute c such that (n/p)^c >= p
			size_t c = 1;
			for( ; c < p; c *= np ) {;}
			//store splitting stage
			k1 = n / c;
		}

		/**
		 * Default constructor. Initialises with meaningless data.
		 * Only to be used internally as a placeholder.
		 */
		MultiBSP_FFT_data() : n( 0 ), k1( 0 ), superstep( 0 ), sign( 0 ), MultiBSP_FFT_data_common() {}

		/**
		 * Direct constructor. For internal use only.
		 *
		 * @param _n    The local problem size.
		 * @param _k1   The splitting stage.
		 * @param _sign The mode of the FFT.
		 * @param _s    The PID corresponding to this internal node.
		 * @param _leaf The parent's leaf ID.
		 * @param _PIDs The PIDs of the parent node.
		 * @param _Ps   The number of processes on each parent level.
		 */
		MultiBSP_FFT_data( const size_t _n, const size_t _k1, const signed char _sign, const size_t _s,
					const size_t _leaf, const std::vector< size_t > _PIDs,
					const std::vector< size_t > _Ps ) :
			n( _n ), k1( _k1 ), superstep( 0 ), sign( _sign ),
			PIDs( _PIDs ), Ps( _Ps ), MultiBSP_FFT_data_common() {
			//append current ID
			PIDs.push_back( _s );
			//append current _P
			Ps.push_back( _P );
			//derive leaf ID
			leaf_id = _leaf + _s * MultiBSP_Computer< _P, _M, _tail... >::processors();
			//initialise child structures
			for( size_t s = 0; s < _P; ++s ) {
				//note we should pass parameters valid for the child node
				children[ s ] = MultiBSP_FFT_internal< _P, _M, _tail... >( n / _P, k1, sign, s, leaf_id, PIDs, Ps, this );
			}
		}

	public:

		/** Local vector length. */
		size_t n;

		/** Splitting point. */
		size_t k1;

		/** Current superstep. */
		size_t superstep;

		/** Leaf-level ID. Assumes a 0-ID on all subtree levels. */
		size_t leaf_id;

		/** Mode of the FFT. */
		signed char sign;

		/** Process IDs. */
		std::vector< size_t > PIDs;

		/** Number of processes on each parent level. */
		std::vector< size_t > Ps;

		/**
		 * Base constructor. Automatically derives k1.
		 *
		 * @param _n    The local problem size.
		 * @param _sign The mode of the FFT.
		 */
		MultiBSP_FFT_data( const size_t _n, const signed char _sign ) :
			n( _n ), superstep( 0 ), leaf_id( 0 ), sign( _sign ), MultiBSP_FFT_data_common()  {
			//initialise k1
			k1_init();
			//append _P
			Ps.push_back( _P );
			//initialise child structures
			for( size_t s = 0; s < _P; ++s ) {
				//note we should pass parameters valid for the child node
				children[ s ] = MultiBSP_FFT_internal< _P, _M, _tail... >( n / _P, k1, sign, s, leaf_id, PIDs, Ps, this );
			}
		}

		/**
		 * Initialises this internal-level FFT data.
		 * Should be called from within an SPMD section.
		 */
		virtual void initialise() {
			//initialise P-bitreversal; only if we are top-level
			if( superstep == 0 ) {
				//derive P
				const size_t P = bsp_nprocs() * MultiBSP_Computer< _P, _M, _tail... >::processors();
				//allocate array
				bitrev_p = new size_t[ P ];
				//initialise values for bitreversal
				bitrev_init( P, bitrev_p );
			}
		}

		/**
		 * Destroys this internal-level FFT data.
		 */
		virtual void destroy() {
			//guard against empty array
			if( bitrev_p != NULL ) {
				//free
				delete [] bitrev_p;
				//set to NULL
				bitrev_p = NULL;
			}
		}

		/** @see Initialisable::retrieve( id ) */
		virtual MultiBSP_FFT_internal< _P, _M, _tail... > & retrieve( const size_t id = bsp_pid() ) {
			return children[ id ];
		}
};

/**
 * This is the MultiBSP_FFT implementation on leaf nodes.
 *
 * @param _P The total number of processes on this level of the MultiBSP computer.
 * @param _M The amount of memory available on this level of the computer.
 * @param _tail The remainder MultiBSP computer paramters (assumed empty here).
 */
template< size_t _P, size_t _M, size_t ... _tail >
class MultiBSP_FFT : public mcbsp::MultiBSP_program< _P, _M, _tail... > {

	protected:

		/** Local data. */
		MultiBSP_FFT_internal< _P, _M, _tail... > &data;

		/** Input vector. */
		Distributed_Vector< _P, _M, double, _tail... > &x;

		/** Whether our process ID is in bit-reverse. */
		bool bitrev;

		/**
		 * Code for swap-based permutations.
		 * Warning: if the permutation array is not reducible to n
		 * swaps, the permutation will fail!
		 *
		 * @param x The array which to permute.
		 * @param n The size of the array.
		 * @param s The permutation array of length n.
		 */
		static void permute( double * __restrict__ const x, const size_t n, const size_t * __restrict__ const s ) {
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
		void redistribute() {
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
			assert( data.n % _P == 0 );
			//get my process ID, account for possible reversal
			//FIXME TODO const size_t s = data.superstep == 0 ? data.bitrev_p[ bsp_pid() ] : bsp_pid();
			const size_t s = bsp_pid();
			//get a raw handle of x
			double * __restrict__ const raw_x = &(x[0]);
			//derive the size of each message to send to the other processes (in #doubles)
			const size_t local_size = 2 * data.n / _P;
			//allocate buffer
			double * __restrict__ const buffer = data.buffer;
			//initialise TODO move to persistent storage
			bsp_push_reg( raw_x, 2*data.n );
			bsp_sync();
			//index in the buffer array
			size_t j = 0;
			//for each other process
			for( size_t k = 0; k < _P; ++k ) {
				//remember chunk start
				double * const chunk_start = buffer + j;
				//build the message by following the
				//cyclic distribution locally
				for( size_t i = k; i < data.n; i += _P ) {
					buffer[ j++ ] = raw_x[   2*i   ];
					buffer[ j++ ] = raw_x[ 2*i + 1 ]; 
				}
				//calculate offset at process k
				const size_t offset = data.superstep == 0 ? data.bitrev_p[ s ] * local_size : s * local_size;
				//send chunk to process k using bsp_hpput
				//for efficient buffering.
				bsp_hpput( k, chunk_start, raw_x, offset, local_size );
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

		/** Apply set of twiddles to the local x. */
		void twiddle() {
			//retrieve handle to raw x data
			double * __restrict__ const raw_x = &(x[0]);

			//sanity check
			assert( data.sign != 0 );

			//switch signs (forward if > 0, backward if < 0)
			if( data.sign > 0 ) {
				//loop over all local elements
				for( size_t j = 0; j < data.n; ++j ) {
					//get aliae of compute elements
					const double real_weight = data.tw[   2*j   ];
					const double imag_weight = data.tw[ 2*j + 1 ];
					const double real_x      =   raw_x[   2*j   ];
					const double imag_x      =   raw_x[ 2*j + 1 ];
					//write back twiddled factors
					x[   2*j   ] = real_weight * real_x - imag_weight * imag_x;
					x[ 2*j + 1 ] = imag_weight * real_x + real_weight * imag_x;
				}
			} else {
				//data.mode > 0
				//loop over all local elements
				for( size_t j = 0; j < data.n; ++j ) {
					//get aliae of compute elements
					const double real_weight =  data.tw[   2*j   ];
					//(account for sign change)
					const double imag_weight = -data.tw[ 2*j + 1 ];
					const double real_x      =    raw_x[   2*j   ];
					const double imag_x      =    raw_x[ 2*j + 1 ];
					//write back twiddled factors
					x[   2*j   ] = real_weight * real_x - imag_weight * imag_x;
					x[ 2*j + 1 ] = imag_weight * real_x + real_weight * imag_x;
				}
			}
			//done
		}

		//DBG
		void printVec( const char * const str ) {
			std::cout << "*** " << data.leaf_id << "=( " << data.PIDs[ 0 ];
			for( size_t k = 1; k < data.PIDs.size(); ++k ) {
				std::cout << ", " << data.PIDs[ k ];
			}
			std::cout << " ); " << str;
			for( size_t i = 0; i < 2 * data.n; i += 2 ) {
				std::cout << x[ i ] << "+" << x[i+1] << "i ";
			}
			std::cout << std::endl;
		}

		/** Generate input data */
		void generateInput() {
			size_t cumulPs[ data.Ps.size() ];
			//get number of processes on each level
			cumulPs[ data.PIDs.size() - 1 ] = 1;
			for( size_t i = data.PIDs.size()-2; i < data.PIDs.size(); --i ) {
				cumulPs[ i ] = data.Ps[ i+1 ] * cumulPs[ i+1 ];
			}
			//get flat process index
			size_t flatP = 0;
			for( size_t i = 0; i < data.PIDs.size(); ++i ) {
				flatP += data.PIDs[ i ] * cumulPs[ i ];
			}
			//initialise x
			for( size_t i = 0; i < data.n; ++i  ) {
				//get global position
				const size_t index = i * cumulPs[ 0 ] * data.Ps[ 0 ] + flatP;
				x[   2*i   ] = index;
				x[ 2*i + 1 ] = 0;
			}
		}

		/** Stage 1 code. */
		void stage1() {
			//DBG
			//printVec( "Pre-stage1 local vector: " );

			//sanity check
			assert( 2 * data.n == x.length() );
			//first do bit-reversal, followed by an undo last stages of the bit-reversion;
			//FFTW does not have unordered FFTs, so a call to FFTW would re-do the last
			//stages. We don't want that. First check if this entire sequence is necessary:
			if( data.k1 < data.n ) { //otherwise, k1 == n and no action is required whatsoever.
				//bit-reverse input
				permute( &(x[0]), data.n, data.bitrev_n );
				//undo last stages of the bit-reversion
				for( size_t r = 0; r < data.n / data.k1 && data.k1 > 2; ++r ) { //if k1<2, the bit-reversal is an identity operation
					//calculate offset of this chunk
					const size_t offset = 2 * r * data.k1;
					permute( &(x[offset]), data.k1, data.bitrev_k1 );
				}
			} else {
				//sanity check
				assert( data.k1 == data.n );
			}

			//get a pointer to the raw array representation of x
			double * __restrict__ const raw_x = &(x[0]);
			//delegate n/k1 FFTs of length k1 to FFTW
			fftw_execute_dft( data.fftw_stage1, (fftw_complex*)(raw_x), (fftw_complex*)(raw_x) );

			//DBG
			//printVec( "Post-stage1 local vector: " );

			//done with leaf-level compute stage 1; first proceed with global redistribution.
		}

		/** Stage 2 code. */
		void stage2() {
			//do redistribution from block to cyclic
			redistribute();

			//get a raw handle of x
			double * __restrict__ raw_x = &(x[0]);

			//do a bit-reversion on x to reduce UGFFT to GFFT
			//TODO: can this be integrated into redistribute?
			permute( raw_x, data.n, data.bitrev_n );

			//DBG
			//printVec( "Pre-twiddle local vector: " );

			//apply twiddle to reduce GFFT to FFT
			twiddle();

			//DBG
			//printVec( "After-twiddle local vector: " );

			//apply final stage2 FFT
			fftw_execute_dft( data.fftw_stage2, (fftw_complex*)(raw_x), (fftw_complex*)(raw_x) );

			//DBG
			//printVec( "After stage-2 local vector: " );
		}

		/** The leaf-case parallel code. */
		virtual void spmd() {
			//if we control the execution flow
			if( data.superstep == 0 ) {
				//generate input vector
				generateInput();
				//prepare benchmark
				double start, end;
				double times[ bsp_nprocs() * 10 ];
				bsp_push_reg( times, bsp_nprocs() * 10 );
				bsp_sync();
				//execute benchmark
				for( size_t rep1 = 0; rep1 < 10; ++rep1 ) {
					//get start time
					start = bsp_time();
					for( size_t rep2 = 0; rep2 < 30; ++rep2 ) {
						//do stage 1
						stage1();
						//do stage 2
						stage2();
					}
					//get end time
					end = bsp_time();
					//store
					times[ bsp_pid() * 10 + rep1 ] = (end - start) / static_cast< double >( 30 );
				}
				//communicate time taken
				for( size_t k = 0; k < bsp_nprocs(); ++k )
					if( k == bsp_pid() )
						bsp_put( k, times + bsp_pid() * 10, times, bsp_pid() * 10, 10 );
				bsp_sync();
				//derive statistics
				double average = 0.0;
				double minimum = INFINITY;
				for( size_t rep1 = 0; rep1 < 10; ++rep1 ) {
					//derive slowest times over all processes
					for( size_t q = 1; q < bsp_nprocs(); ++q ) {
						const double cur = times[ q * 10 + rep1 ];
						if( cur > times[ rep1 ] )
							times[ rep1 ] = cur;
					}
					//calculate average, minimum on the fly as well
					average += times[ rep1 ] / static_cast< double >( 10 );
					if( times[ rep1 ] < minimum )
						minimum = times[ rep1 ];
				}
				//calculate variance
				double variance = 0.0;
				for( size_t rep1 = 0; rep1 < 10; ++rep1 ) {
					const double diff = average - times[ rep1 ];
					variance += diff * diff / static_cast< double >( 10 - 1 );
				}
				//calculate stddev, flop count
				const double n = static_cast< double >( data.n * _P );
				const double stddev = sqrt( variance );
				const double nflops = 5*n*log(n) / log( 2.0 ) + 2*n;
				const double Gflops = nflops / 1000000000.0;
				//output results
				if( bsp_pid() == 0 ) {
					printf( "Average time per FFT: %lf ms., minimum time measured: %lf ms., stddev: %lf ms.\n", (average*1000.0), (minimum*1000.0), (stddev*1000.0) );
					printf( "Average computing rate in FFT: %lf Gflop/s, maximum flop rate: %lf Gflop/s.\n", (Gflops/average), (Gflops/minimum) );
				}
				//done
				return;
			}
			//otherwise, switch execution path according to current stage
			if( data.superstep == 1 ) {
				generateInput();
			} else if( data.superstep % 2 == 0 ) {
				stage1();
			} else if( data.superstep % 2 == 1 ) {
				stage2();
			} else {
				bsp_abort( "Arrived at undefined superstep %ld!\n", data.superstep );
			}
			//increment current superstep
			++(data.superstep);
		}

		/**
		 * Function for supporting mcbsp::MultiBSP_program::begin().
		 *
		 * @return A new instance of this class for a sibling process to work on.
		 */
		virtual MultiBSP_FFT< _P, _M, _tail... > * newInstance() {
			return new MultiBSP_FFT< _P, _M, _tail... >( data.parent().retrieve(), x.parent().retrieve() );
		}

	public:

		/**
		 * Direct constructor. For internal use only.
		 *
		 * @param _data Local FFT data.
		 */
		MultiBSP_FFT( MultiBSP_FFT_internal< _P, _M, _tail... > &_data, Distributed_Vector< _P, _M, double, _tail... > &_x ) : data( _data ), x( _x ) {}

		/**
		 * Base constructor.
		 *
		 * @param _data Global FFT data.
		 */
		MultiBSP_FFT( MultiBSP_FFT_data< _P, _M, _tail... > &_data, Distributed_Vector_Root< double, _P, _M, _tail... > &_x ) : data( _data.retrieve( 0 ) ), x( _x.retrieve( 0 ) ) {}
};

template< size_t _P, size_t _M, size_t _subP, size_t _subM, size_t ... _tail >
class MultiBSP_FFT< _P, _M, _subP, _subM, _tail... > : public mcbsp::MultiBSP_program< _P, _M, _subP, _subM, _tail... > {

	protected:

		//local data
		MultiBSP_FFT_internal< _P, _M, _subP, _subM, _tail... > &data;

		//input vector to work on
		Distributed_Vector< _P, _M, double, _subP, _subM, _tail... > &x;

		/** Performs top-level communication for redistribution. */
		void global_redistribute() {
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

			//our ID
			const size_t s = bsp_pid();
			//number of processors in this subtree
			const size_t subtree_P = MultiBSP_Computer< _subP, _subM, _tail... >::processors();
			//get leaf-level size
			const size_t np = data.n / subtree_P;
			//array of destination indices TODO move to persistent data
			size_t * __restrict__ const dest_ind = new size_t[ np ];
			//index counter
			size_t index = 0;
			//loop over blocks of size subtree_P
			for( size_t i = 0; i < data.n; i += subtree_P ) {
				//get lower bits
				const size_t lower = i % np;
				//get upper bits; if we are the top-level run,
				//account for bit-reversed process IDs here.
				//Then redistribution automagically fixes
				//the process IDs.
				const size_t upper_bits = i / np + s * subtree_P;
				const size_t upper = data.superstep == 0 ? data.bitrev_p[ upper_bits ] : upper_bits;
				//derive global index
				const size_t global = upper * np + lower;
				//derive destination process
				const size_t destproc = (global / subtree_P) % _P;
				//derive destination index
				dest_ind[ index ] = ((global / subtree_P / _P) * subtree_P +
							(global % subtree_P)) % data.n;
				//check if destination is local
				if( destproc == s && i == dest_ind[ index ] ) {
					//no action required
					continue;
				}
				//get pointer to local block
				const double * __restrict__ const block = x.pointerToElement( 2 * i );
				//send message
				bsp_hpsend( destproc, dest_ind + index, block, 2 * subtree_P );
				//increment index
				++index;
			}
			//execute communication
			bsp_sync();
			//define pointers to tag, payloads
			size_t * tag;
			double * payload;
			//handle communication
			while( bsp_hpmove( (void **)(&tag), (void**)(&payload) ) != SIZE_MAX ) {
				//get pointer to block to be overwritten
				double * __restrict__ const block = x.pointerToElement( 2 * *tag );
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
			/*for( size_t k = 0; k < _P; ++k ) {
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

		//the parallel code at this level
		virtual void spmd() {
			//initialise (TODO: move this to persistent data
			bsp_set_tagsize< size_t >();
			bsp_sync();

			//check if we have the control flow
			if( data.superstep == 0 ) {
				//do initialisation of input
				this->bsp_recurse();
				//prepare benchmark
				double start, end;
				double times[ bsp_nprocs() * 10 ];
				bsp_push_reg( times, bsp_nprocs() * 10 );
				bsp_sync();
				//execute benchmark
				for( size_t rep1 = 0; rep1 < 10; ++rep1 ) {
					//get start time
					start = bsp_time();
					for( size_t rep2 = 0; rep2 < 30; ++rep2 ) {
						//do stage1
						this->bsp_recurse();
						//do redistribute
						global_redistribute();
						//do lower-level redistributes and stage2s
						this->bsp_recurse();
					}
					//get end time
					end = bsp_time();
					//store
					times[ bsp_pid() * 10 + rep1 ] = (end - start) / static_cast< double >( 30 );
				}
				//communicate time taken
				for( size_t k = 0; k < bsp_nprocs(); ++k )
					if( k == bsp_pid() )
						bsp_put( k, times + bsp_pid() * 10, times, bsp_pid() * 10, 10 );
				bsp_sync();
				//derive statistics
				double average = 0.0;
				double minimum = INFINITY;
				for( size_t rep1 = 0; rep1 < 10; ++rep1 ) {
					//derive slowest times over all processes
					for( size_t q = 1; q < bsp_nprocs(); ++q ) {
						const double cur = times[ q * 10 + rep1 ];
						if( cur > times[ rep1 ] )
							times[ rep1 ] = cur;
					}
					//calculate average, minimum on the fly as well
					average += times[ rep1 ] / static_cast< double >( 10 );
					if( times[ rep1 ] < minimum )
						minimum = times[ rep1 ];
				}
				//calculate variance
				double variance = 0.0;
				for( size_t rep1 = 0; rep1 < 10; ++rep1 ) {
					const double diff = average - times[ rep1 ];
					variance += diff * diff / static_cast< double >( 10 - 1 );
				}
				//calculate stddev, flop count
				const double n = static_cast< double >( data.n * _P );
				const double stddev = sqrt( variance );
				const double nflops = 5*n*log(n) / log( 2.0 ) + 2*n;
				const double Gflops = nflops / 1000000000.0;
				//output results
				if( bsp_pid() == 0 ) {
					printf( "Average time per FFT: %lf ms., minimum time measured: %lf ms., stddev: %lf ms.\n", (average*1000.0), (minimum*1000.0), (stddev*1000.0) );
					printf( "Average computing rate in FFT: %lf Gflop/s, maximum flop rate: %lf Gflop/s.\n", (Gflops/average), (Gflops/minimum) );
				}
				//done
				return;
			}
			//adapt behaviour according to current superstep
			if( data.superstep == 1 ) {
				//just recurse; input init happens at leafs
				this->bsp_recurse();
			} else if( data.superstep % 2 == 0 ) {
				//just recurse; no global actions
				this->bsp_recurse();
			} else if( data.superstep % 2 == 1 ) {
				//first do global redistribute
				global_redistribute();
				//then recurse for completion
				this->bsp_recurse();
			}
			//increase current superstep
			++(data.superstep);
		}

		virtual MultiBSP_FFT< _P, _M, _subP, _subM, _tail... > * newInstance() {
			return new MultiBSP_FFT< _P, _M, _subP, _subM, _tail... >( data.parent().retrieve(), x.parent().retrieve() );
		}

		virtual MultiBSP_FFT< _subP, _subM, _tail... > * newChild( const size_t s = bsp_pid() ) {
			//if this is not a repeated run
			if( data.retrieve( s ).superstep == 0 ) {
				//signal it should start at the first superstep
				data.retrieve( s ).superstep = 1;
			}
			//return new child
			return new MultiBSP_FFT< _subP, _subM, _tail... >( data.retrieve( s ), x.retrieve( s ) );
		}

	public:

		/**
		 * Direct constructor. For internal use only.
		 *
		 * @param _data Local FFT data.
		 */
		MultiBSP_FFT( MultiBSP_FFT_internal< _P, _M, _subP, _subM, _tail... > &_data, Distributed_Vector< _P, _M, double, _subP, _subM, _tail... > &_x ) : data( _data ), x( _x ) {}

		/**
		 * Base constructor.
		 *
		 * @param _data Global FFT data.
		 */
		MultiBSP_FFT( MultiBSP_FFT_data< _P, _M, _subP, _subM, _tail... > &_data, Distributed_Vector_Root< double, _P, _M, _subP, _subM, _tail... > &_x ) : data( _data.retrieve( 0 ) ), x( _x.retrieve( 0 ) ) {}
};

