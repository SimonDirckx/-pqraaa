
#ifndef _H_MCBSP_ALLREDUCE
#define _H_MCBSP_ALLREDUCE

#include <limits>

#include "mcbsp-templates.hpp"
#include "multibsp.hpp"
#include "templates/collectives/reduce.hpp"
#include "templates/collectives/broadcast.hpp"

namespace mcbsp {
namespace collectives {

	/**
	 * Codes an allreduce collective operation.
	 * TODO further documentation goes here
	 */
	template< typename Operation, typename ValueType >
	class AllReduce :
		public virtual MultiBSP_template,
		protected virtual Reduce< Operation, ValueType >,
		protected virtual Broadcast< ValueType > {
	
		template< typename T >
		friend class Dot;

		private:

		protected:

			void init() {
				Broadcast< ValueType >::init();
				Reduce< Operation, ValueType >::init();
			}

			void destroy() {
				Broadcast< ValueType >::destroy();
				Reduce< Operation, ValueType >::destroy();
			}

			ValueType leaf( const ValueType * const array = NULL, const size_t size = 0 ) {
				//set values
				Reduce< Operation, ValueType >::to_reduce = array;
				Reduce< Operation, ValueType >::length    = size;
				Reduce< Operation, ValueType >::target    = false;
				Broadcast< ValueType >::size = 1;
				//do local reduce
				Reduce< Operation, ValueType >::localReduce();
				//do communicate reduced values
				Reduce< Operation, ValueType >::all_to_one_reduce();
				//let ancestors do upper-level reductions
				bsp_up();
				//pass reduced result into broadcast template
				Broadcast< ValueType >::receive_buffer[ 0 ] = Reduce< Operation, ValueType >::reduced;
				//broadcast global result
				Broadcast< ValueType >::bcast( &(Reduce< Operation, ValueType >::reduced) );
				//allow communication to finish
				bsp_sync();
				//return broadcasted value
				return *(Broadcast< ValueType >::receive_buffer);
			}

			void nested() {
				//communicate
				Reduce< Operation, ValueType >::all_to_one_reduce();
				//let ancestors do upper-level reductions
				bsp_up();
				//pass reduced result into broadcast template
				Broadcast< ValueType >::receive_buffer[ 0 ] = Reduce< Operation, ValueType >::reduced;
				//broadcast global result
				Broadcast< ValueType >::bcast( &( Reduce< Operation, ValueType >::reduced) );
				//go back down
				bsp_down();
			}

			ValueType apply( const ValueType * const array = NULL , const size_t size = 0 ) {
				//dispatch
				if( bsp_leaf() ) {
					return leaf( array, size );
				} else {
					//nested branch
					nested();
					//return neutral value in nested case
					return Operation::zero;
				}
			}

		public:

	};

} //collectives
} //mcbsp

#endif

