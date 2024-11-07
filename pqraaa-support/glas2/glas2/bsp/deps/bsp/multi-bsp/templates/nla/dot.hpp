
#ifndef _H_MCBSP_DOT
#define _H_MCBSP_DOT

#include <assert.h>

#include "ops/ops.hpp"
#include "multibsp.hpp"
#include "mcbsp-templates.hpp"
#include "collectives/allreduce.hpp"

namespace mcbsp {
namespace nla {

	/**
	 * Template which computes the dot product (inner product) of two arrays.
	 *
	 * By default, uses the standard semi-ring (+,*). This can be adapted by
	 * modifying the template parameters `Addition' and `Multiplication'.
	 * Makes use of the Allreduce collective.
	 *
	 * @tparam ValueType      The type of entries to reduce.
	 * @tparam Addition       The addition action of the semi-ring under which
	 * 			  to compute this inner product. The action must
	 * 			  operate on values of type ValueType.
	 * @tparam Multiplication The multiplicative action of the semi-ring under
	 * 			  which to compute this inner product. The action
	 * 			  must operate on values of type ValueType.
	 */
	template< typename ValueType, class Addition = mcbsp::ops::Add< ValueType >, class Multiplication = mcbsp::ops::Mul< ValueType > >
	class Dot :
		public virtual MultiBSP_template,
		protected virtual mcbsp::collectives::AllReduce< mcbsp::ops::Add< ValueType >, ValueType >
	{

		private:

			/** Convenience typedef for the AllReduce template used by this template. */
			typedef mcbsp::collectives::AllReduce< Addition, ValueType > AllReduce;

			/** The zero element corresponding to the standard multiplication on the current semi-ring. */
			static constexpr ValueType zero = Addition::identity;

		protected:

			/** Initialises this template for use. */
			void init() {
				AllReduce::init();
			}

			/** Destroys this template. */
			void destroy() {
				AllReduce::destroy();
			}

			/**
			 * Leaf SPMD section code.
			 *
			 * @param x Local view of one of the vectors to take a dot product with.
			 * @param y Local view of one of the vectors to take a dot product with.
			 * @param n The length of both x and y.
			 * @return The inner product \f$ \alpha = x^Ty = \sum_{i=0}^n x_iy_i \f$.
			 */
			ValueType leaf( const ValueType *__restrict__ const x, const ValueType *__restrict__ const y, const size_t n ) {
				//calculate local contribution; set initial value, identity under addition (i.e., zero):
				ValueType alpha = Addition::identity;
				//sanity check
				assert( n == 0 || (x != NULL && y != NULL) );
				//add contributions of both vectors element-by-element
				for( size_t i = 0; i < n; ++i ) {
					//The below performs `alpha += x[ i ] * y[ i ];' in a generic way. First,
					//copy one vector element to apply the multiplication on:
					ValueType val = x[i];
					//then apply multiplication action
					Multiplication::apply( val, y[i] );
					//and finally apply the addition action
					Addition::apply( alpha, val );
				}
				//now do allreduce on local contribution to the inner product,
				//and return the resulting global inner product.
				return AllReduce::leaf( &alpha, 1 );
			} 

			/**
			 * The nested SPMD section code.
			 */
			void nested() {
				//in the internal node code, simply do allreduce
				AllReduce::nested();
			}

			/**
			 * The generic entry point; decides automatically whether to dispatch to
			 * the leaf or nested codes. It is slightly more efficient to directly
			 * call the correct SPMD section. Options are like that of the leaf
			 * function call, but are optional.
			 *
			 * @see leaf.
			 * @return On a leaf node, return the global inner product;
			 *         on a nested node, returns 0.
			 */
			ValueType apply( const ValueType *__restrict__ const x = NULL, const ValueType *__restrict__ const y = NULL, const size_t n = 0 ) {
				//leaf/nested case branch
				if( bsp_leaf() ) {
					//leaf dispatch
					return Dot< ValueType >::leaf( x, y, n );
				} else {
					//nested dispatch
					Dot< ValueType >::nested();
					//nested case returns zero
					return zero;
				}
				//done
			}

	};

} //nla
} //mcbsp

#endif

