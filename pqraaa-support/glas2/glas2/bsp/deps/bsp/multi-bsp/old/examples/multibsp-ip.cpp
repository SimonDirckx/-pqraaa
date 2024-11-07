
#include "multibsp-ip.hpp"

#define MULTIBSP_TOP_P 1
#define MULTIBSP_COMPUTER MULTIBSP_TOP_P, 1l<<33, 4, 1l<<25, 1, 1l<<15

//alternative computer model
//#define MULTIBSP_TOP_P 2
//#define MULTIBSP_COMPUTER MULTIBSP_TOP_P, 1l<<33, 2, 1l<<25, 1, 1l<<15

int main() {

	const size_t N = 1l<<18;

	std::cout << "Compiled with MultiBSP computer model:\n";
	MultiBSP_Computer< MULTIBSP_COMPUTER >::print();

	Distributed_Vector_Root< double, MULTIBSP_COMPUTER > x( N );
	Distributed_Vector_Root< double, MULTIBSP_COMPUTER > y( N );
	Replicated_Vector_Root< double, MULTIBSP_COMPUTER > alpha( 1 );

	std::vector< Initialisable< MULTIBSP_COMPUTER > * > variables;
	variables.push_back( &x );
	variables.push_back( &y );
	variables.push_back( &alpha );

	Initialiser< MULTIBSP_COMPUTER > initialiser = Initialiser< MULTIBSP_COMPUTER >( variables );
	initialiser.begin( MULTIBSP_TOP_P );
	
	MultiBSP_IP< MULTIBSP_COMPUTER > ip( x, y, alpha );
	ip.begin();

	std::cout << "Result of inner product calculation: " << alpha.retrieve( 0 )[ 0 ] << std::endl;
	std::cout << "              Expectation value was: " << (N/4.0) << std::endl;

	Destructor< MULTIBSP_COMPUTER > destructor = Destructor< MULTIBSP_COMPUTER >( variables );
	destructor.begin( MULTIBSP_TOP_P );

	return 0;
};

