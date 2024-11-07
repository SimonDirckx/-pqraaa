
#include <cstdlib>

void hardc_add( double * left, const double * right, const size_t n ) {
	for( size_t i = 0; i < n; ++i ) {
		left[ i ] += right[ i ];
	}
}

