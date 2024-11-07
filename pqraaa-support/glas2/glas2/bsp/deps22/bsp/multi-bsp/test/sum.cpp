
#include <cstdlib>

#include "blas/blas.hpp"
#include "collectives/broadcast.hpp"

#ifndef MULTIBSP_COMPUTER
 #define MULTIBSP_COMPUTER 2, 1l<<33, 2, 1l<<25, 1, 1l<<15
#endif

/** External function doing an optimised vector addition. */
extern void hardc_add( double * __restrict__ x, const double * __restrict__ y, const size_t n );

/** External function doing a (hopefully) optimised vector addition using templated code. */
extern void templ_add( double * __restrict__ x, const double * __restrict__ y, const size_t n );

/** Test Multi-BSP algorithm which sums two vectors using various methods. */
class Test :
	public virtual mcbsp::MultiBSP_algorithm< Test >,
	protected virtual mcbsp::collectives::Broadcast< size_t >
{

	private:

		/** Convenience typedef to reach functions corresponding to the Broadcast template. */
		typedef mcbsp::collectives::Broadcast< size_t > Broadcast_SIZE_T;

		/** One of the two input vectors. */
		double *x;

		/** One of the two output vectors. */
		double *y;

		/** Leaf SPMD code. */
		virtual void leaf_spmd() {
			//initialise templates
			Broadcast_SIZE_T::init();
			bsp_up();

			//algorithm initialisation: choose N
			size_t N;
			if( bsp_lid() == 0 ) {
				srand( 1337 );
				N = static_cast< size_t >( rand() ) % 100000000ul;
				std::cout << "0: randomly selected N = " << N << "(to completely prevent benchmark cheating by overly smart compilers)" << std::endl;
				Broadcast_SIZE_T::leaf( &N, 1 );
			} else {
				N = *(Broadcast_SIZE_T::leaf( 1 ));
			}

			//allocate, initialise arrays
			std::cout << bsp_lid() << ": will run test using a local array of size " << N << std::endl;
			x = new double[ N ];
			y = new double[ N ];
			for( size_t i = 0; i < N; ++i ) {
				x[ i ] = y[ i ] = 1.0;
			}

			//verification
			mcbsp::BLAS< double >::axpy( y, x, N );
			for( size_t i = 0; i < N; ++i ) {
				if( y[ i ] != 2.0 ) {
					bsp_abort( "Verification failed at entry %ul (%lf instead of 2.0)\n", i, x[i] );
				}
			}
			//warm-up kernels
			hardc_add( x, y, N );
			templ_add( x, y, N );
			//sync and do benchmarks
			bsp_up();

			//hardcoded addition
			double hardcoded = bsp_time();
			hardc_add( x, y, N );
			hardcoded = bsp_time() - hardcoded;
			std::cout << bsp_lid() << ": hardcoded version finished in " << hardcoded << " seconds." << std::endl;
			bsp_up();

			//partially templated
			double partially = bsp_time();
			templ_add( x, y, N );
			partially = bsp_time() - partially;
			std::cout << bsp_lid() << ": templated version  (ops::Add)  finished in " << partially << " seconds." << std::endl;
			bsp_up();

			//templated addition
			double templated = bsp_time();
			mcbsp::BLAS< double >::axpy( x, y, N );
			templated = bsp_time() - templated;
			std::cout << bsp_lid() << ": templated version (BLAS::axpy) finished in " << templated << " seconds." << std::endl;

			//destroy templates
			Broadcast_SIZE_T::destroy();
			bsp_up();
		}

		/** Nested SPMD code. */
		virtual void nested_spmd() {
			//initialise templates
			Broadcast_SIZE_T::init();
			bsp_up();
			bsp_down();

			//do broadcast of N
			Broadcast_SIZE_T::nested();

			//sync after verification
			bsp_up();
			bsp_down();

			//sync after hardcoded addition benchmark
			bsp_up();
			bsp_down();

			//sync after templated addition benchmark
			bsp_up();
			bsp_down();

			//destroy templates
			Broadcast_SIZE_T::destroy();
			bsp_up();
			bsp_down();
		}

};

/** Test program entry point. */
int main() {
	//report input parameters
	std::cout << "Configured using the following MultiBSP computer:\n";
	std::cout << mcbsp::MultiBSP_computer< MULTIBSP_COMPUTER >::getDescription() << std::endl;
	//do test
	Test test;
	mcbsp::MultiBSP_computer< MULTIBSP_COMPUTER >::bsp_begin( test );
}

