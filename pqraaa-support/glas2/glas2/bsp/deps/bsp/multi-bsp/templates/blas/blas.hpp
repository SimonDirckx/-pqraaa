
#ifndef _H_MCBSP_BLAS
#define _H_MCBSP_BLAS

#include "ops/ops.hpp"
#include "multibsp.hpp"

namespace mcbsp {

/**
 * Defines the BLAS for generic semi-rings. By default, will take the standard
 * ring under addition and multiplication, i.e., (+,*). The type of object
 * operated on, is templated as well.
 *
 * The rationale of this BLAS library is as follows:
 * 	-Maximum genericity,
 * 	-maximum use of default arguments,
 * 	-maximum performance through template specialisation,
 * 	-the base function names are as the de-facto BLAS standard,
 * 	-data modified by BLAS calls are put first in the argument list,
 * 	-parameters with common default arguments are placed last in decreasing common-ness,
 *	-ties in argument placing are resolved by dimensions (matrices > vectors > scalars) first,
 *	-and by alphabet second.
 *
 * Currently implemented BLAS-1 functions:
 * 	-axpy.
 * 
 * Currently implemented BLAS-2 functions:
 * 	-none.
 *
 * Currently implemented BLAS-3 functions:
 * 	-none.
 *
 * @tparam ValueType The type of values which represent values from the set
 *                   the semi-ring is defined on.
 * @tparam Addition  The additive action of the semi-ring.
 * @tparam Multiplication The multiplicative action of the semi-ring.
 *
 * @see mcbsp::ops::Add (the default additive action)
 * @see mcbsp::ops::Mul (the default multiplicative action)
 */
template< typename ValueType, class Addition = mcbsp::ops::Add< ValueType >, class Multiplication = mcbsp::ops::Mul< ValueType > >
struct BLAS {

	/**
	 * Computes \f$ y := y + \alpha x \f$.
	 */
	static void axpy( ValueType * __restrict__ const y,
		const ValueType * __restrict__ const x,
		const size_t n, const ValueType alpha = Addition::identity,
		const size_t stride_x = 1, const size_t stride_y = 1 )
	{
		//sanity checks
		assert( stride_x > 0 );
		assert( stride_y > 0 );
		//check whether stride is 1
		if( stride_x == 1 && stride_y == 1 ) {
			//check whether applying alpha has no effect on x
			if( alpha == Addition::identity ) {
				//degenerate case: apply additive action of x on y
				Addition::apply( y, x, n );
			} else {
				//do kernel
				for( size_t i = 0; i < n; ++i ) {
					//cache element from x
					ValueType val = x[ i ];
					//apply alpha to val
					Multiplication::apply( val, alpha );
					//add \f$ \alpha x_i \f$ to \f$ y_i. \f$
					Addition::apply( y[ i ], val );
				}
			}
		} else {
			//write out strided axpy. Cache pointers:
			const ValueType * __restrict__ x_p = x;
			      ValueType * __restrict__ y_p = y;
			//check whether applying alpha has no effect on x
			if( alpha == Addition::identity ) {
				//do kernel, a=0
				for( size_t i = 0; i < n; ++i, x_p += stride_x, y_p += stride_y ) {
					Addition::apply( *y_p, *x_p );
				}
			} else {
				//do kernel, with alpha
				for( size_t i = 0; i < n; ++i, x_p += stride_x, y_p += stride_y ) {
					//cache element from x
					ValueType val = *x_p;
					//apply alpha to val
					Multiplication::apply( val, alpha );
					//add \f$ \alpha x_{i\mathit{x\_stride} \f$ to \f$ y_{i\mathit{y\_stride}}. \f$
					Addition::apply( *y_p, val );
				}
			}
		}
	}

}; //BLAS
}  //mcbsp

#endif

