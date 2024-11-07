
#include "multibsp.hpp"
#include "numeric/fft.hpp"

/** The MultiBSP computer this examples assumes */
#define MULTIBSP_COMPUTER 2, 1l<<33, 2, 1l<<25, 1, 1l<<15

/** The local vector length at each leaf BSP process */
#define N (1l<<10)

class MultiBSP_FFT : public mcbsp::MultiBSP_algorithm< MultiBSP_FFT >, protected mcbsp::numeric::FFT {

	private:

		/** Convenience typedef. */
		typedef mcbsp::numeric::FFT FFT;

		/** The local vector to be Fourier transformed. */
		double * x;

		/** Generate input data */
		void generateInput() {
			//allocate x
			if( posix_memalign( (void**)&x, 64, 2 * N * sizeof(double) ) != 0 ) {
				bsp_abort( "Could not allocate aligned array of %ld doubles!\n", N );
			}
			//get #leafs and unique leaf ID
			const size_t totalP = bsp_nleafs();
			const size_t leafID = bsp_lid();
			//initialise x
			for( size_t i = 0; i < N; ++i  ) {
				//get global position
				const size_t index = i * totalP + leafID;
				x[   2*i   ] = index;
				x[ 2*i + 1 ] = 0;
			}
			//register our input vector
			bsp_push_reg( x, 2 * N * sizeof(double) );
		}

		virtual void leaf_spmd() {
			//DBG
			std::cout << bsp_lid() << " enters init.\n";
			//initialise FFT code
			FFT::init( N, 1 );
			//DBG
			std::cout << bsp_lid() << " exits init.\n";
			//initialise input vector
			generateInput();
			//allow init to finish
			bsp_up();
			//DBG
			std::cout << "calling FFT::leaf with input vector @ " << x << "\n";
			//execute a single FFT
			FFT::leaf( x );
			//clean up
			FFT::destroy();
			free( x );
			//done
		}

		virtual void nested_spmd() {
			//initialise
			FFT::init( N, 1 );
			//allow init to finish
			bsp_up();
			//go back down to start doing a single FFT
			bsp_down();
			//execute a single FFT
			FFT::nested( x );
			//clean up
			FFT::destroy();
			//done
		}

};

/** The test code. */
int main() {
	//report input parameters
	std::cout << "Configured using the following MultiBSP computer:\n";
	std::cout << mcbsp::MultiBSP_computer< MULTIBSP_COMPUTER >::getDescription() << std::endl;

	//call algorithm, compile, and run it
	MultiBSP_FFT algorithm;
	mcbsp::MultiBSP_computer< MULTIBSP_COMPUTER >::bsp_begin( algorithm );

	//done
	return 0;
}

