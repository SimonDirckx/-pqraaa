
#include "multibsp-fft.hpp"

//#define MULTIBSP_COMPUTER 1, 1l<<33, 4, 1l<<25, 1, 1l<<15

//DL980
//#define MULTIBSP_COMPUTER 1, 1l<<33, 64, 1l<<32
//#define MULTIBSP_COMPUTER 64, 1l<<33
//#define MULTIBSP_COMPUTER 8, 1l<<33, 8, 1l<<32
//#define MULTIBSP_COMPUTER 2, 1l<<33, 32, 1l<<32
//#define MULTIBSP_COMPUTER 2, 1l<<33, 4, 1l<<32, 8, 1l<<31

//DL580
#define MULTIBSP_COMPUTER 1, 1l<<33, 32, 1l<<32
//#define MULTIBSP_COMPUTER 4, 1l<<33, 8, 1l<<32

int main( int argc, char ** argv ) {

	if( argc != 2 ) {
		std::cout << "Usage: " << argv[ 0 ] << " <m>" << std::endl;
		return EXIT_SUCCESS;
	}

	const size_t n = strtoul( argv[ 1 ], NULL, 0 );
	if( n == 0 || n == SIZE_MAX ) {
		std::cerr << "Erroneous value for m ( " << n << " )!" << std::endl;
		return EXIT_FAILURE;
	}
	const size_t N = 1l << n;

	std::cout << "Compiled with MultiBSP computer model:\n";
	MultiBSP_Computer< MULTIBSP_COMPUTER >::print();

	Distributed_Vector_Root< double, MULTIBSP_COMPUTER > x( 2 * N );
	MultiBSP_FFT_data< MULTIBSP_COMPUTER > fft_data( N, 1 );

	std::vector< Initialisable< MULTIBSP_COMPUTER > * > variables;
	variables.push_back( &x );
	variables.push_back( &fft_data );

	Initialiser< MULTIBSP_COMPUTER > initialiser = Initialiser< MULTIBSP_COMPUTER >( variables );
	initialiser.begin();

	MultiBSP_FFT< MULTIBSP_COMPUTER > fft( fft_data, x );
	fft.begin();

	Destructor< MULTIBSP_COMPUTER > destructor = Destructor< MULTIBSP_COMPUTER >( variables );
	destructor.begin();

	return 0;
};

