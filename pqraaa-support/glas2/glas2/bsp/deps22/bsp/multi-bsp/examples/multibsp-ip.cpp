
#include "mcbsp-templates.hpp"
#include "multibsp.hpp"

#define MULTIBSP_COMPUTER 1, 1l<<33, 4, 1l<<25, 1, 1l<<15
#define N (1l<<18)
#define total_P 4

/**
 * Inner product calculation followed by a gather to the MultiBSP computer root.
 * Direct implementation.
 */
class MultiBSP_IP : public mcbsp::MultiBSP_algorithm< MultiBSP_IP > {

	private:

		double *x, *y;

		virtual void leaf_spmd() {
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
			//compute local inner product
			alpha = 0.0;
			for( size_t i = 0; i < N/total_P; ++i ) {
				alpha += x[ i ] * y[ i ];
			}
			//register buffer
			bsp_push_reg( &alpha, 1 );
			//wait for sibling computations
			bsp_sync();
			//root process derives global inner product, through reduction
			for( unsigned int t = 1; s == 0 && t < P; ++t ) {
				//local buffer
				double buffer;
				//otherwise grab data and compute inner product
 				bsp_direct_get( t, &alpha, 0, &buffer, 1 );
				//reduce
				alpha += buffer;
			}
			//clean up
			delete [] x;
			delete [] y;
			//done
		}

		virtual void nested_spmd() {
			//Our process ID, and number of processes
			const unsigned int s = bsp_pid();
			const unsigned int P = bsp_nprocs();
			//prepare nested run for communication
			bsp_push_reg( &alpha );
			//wait for sibling computations, make valid push_reg
			bsp_sync();
			//root process reduces remote contributions
			for( unsigned int k = 1; s == 0 && k < P; ++k ) {
				//declare buffer
				double buffer;
				//otherwise grab data and compute inner product
				bsp_direct_get( k, &(alpha), 0, &buffer );
				//reduce
				alpha += buffer;
			}
			//done
		}

	public:

		/** Result field. */
		double alpha;

};

int main() {
	//report compilation parameters
	std::cout << "Configured using the following MultiBSP computer:\n";
	std::cout << mcbsp::MultiBSP_computer< MULTIBSP_COMPUTER >::getDescription() << std::endl;

	//call algorithm, compile, and run it
	MultiBSP_IP algorithm;
	mcbsp::MultiBSP_computer< MULTIBSP_COMPUTER >::bsp_begin( algorithm );

	//report statistics
	std::cout << "\nResult of inner product calculation: " << algorithm.alpha << std::endl;
	std::cout <<   "              Expectation value was: " << (static_cast<double>(N)/4.0) << std::endl;

	//done
	return 0;
}

