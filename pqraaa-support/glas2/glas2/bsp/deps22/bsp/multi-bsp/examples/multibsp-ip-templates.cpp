
#include "multibsp.hpp"
#include "nla/dot.hpp"
#include "ops/ops.hpp"

#define MULTIBSP_COMPUTER 2, 1l<<33, 2, 1l<<25, 1, 1l<<15
#define N (1l<<18)
#define total_P 4

/**
 * Inner product calculation followed by a gather to the MultiBSP computer root.
 * Direct implementation.
 */
class MultiBSP_IP : public mcbsp::MultiBSP_algorithm< MultiBSP_IP >, protected mcbsp::nla::Dot< double > {

	private:

		/** Short-hand notation of the Reduce collective operation. */
		typedef mcbsp::nla::Dot< double > Dot;

		/** The vectors of which to compute the inner product. */
		double *x, *y;

		/** The leaf node program. */
		virtual void leaf_spmd() {
			//initialise
			Dot::init();
			//allow initialisation to finish
			bsp_sync();
			//Our process ID, and number of processes
			const unsigned int s = bsp_pid();
			const unsigned int P = bsp_nprocs();
			//create local data
			x = new double[ N/total_P ];
			y = new double[ N/total_P ];
			//initialise local data
			for( size_t i = 0; i < N/total_P; ++i ) {
				x[ i ] = ((double)rand()) / ((double)RAND_MAX);
				y[ i ] = ((double)rand()) / ((double)RAND_MAX);
			}
			//calculate dot product
			alpha = Dot::apply( x, y, N/total_P );
			//clean up
			Dot::destroy();
			delete [] x;
			delete [] y;
			//done
		}

		/** The internal node program. */
		virtual void nested_spmd() {
			//initialise
			Dot::init();
			//allow initialisation to finish
			bsp_sync();
			//Our process ID, and number of processes
			const unsigned int s = bsp_pid();
			const unsigned int P = bsp_nprocs();
			//handle reduction
			Dot::apply();
			//clean up
			Dot::destroy();
			//done
		}

	public:

		/** The calculated inner product. */
		double alpha;

};

/** The test code. */
int main() {
	//report input parameters
	std::cout << "Configured using the following MultiBSP computer:\n";
	std::cout << mcbsp::MultiBSP_computer< MULTIBSP_COMPUTER >::getDescription() << std::endl;

	//call algorithm, compile, and run it
	MultiBSP_IP algorithm;
	mcbsp::MultiBSP_computer< MULTIBSP_COMPUTER >::bsp_begin( algorithm );

	//report statistics
	std::cout << "\nResult of inner product calculation: " << algorithm.MultiBSP_IP::alpha << std::endl;
	std::cout <<   "              Expectation value was: " << (static_cast<double>(N)/4.0) << std::endl;

	//done
	return 0;
}

