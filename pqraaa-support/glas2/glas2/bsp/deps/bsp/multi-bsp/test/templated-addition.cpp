
#include <cstdlib>

#include "ops/ops.hpp"

void templ_add( double * __restrict__ left, const double * __restrict__ right, const size_t n ) {
	mcbsp::ops::Add< double >::apply( left, right, n );
}

