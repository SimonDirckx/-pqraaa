
#include <limits>

#include "mcbsp-templates.hpp"
#include "multibsp.hpp"

#ifndef _H_MCBSP_REDUCE
#define _H_MCBSP_REDUCE

namespace mcbsp {
namespace collectives {

	/**
	 * Codes a reduce collective operation. The final reduction will
	 * be held in a root leaf node, as specifically indicated by the
	 * user.
	 *
	 * @tparam Operation The operator used in reduction.
	 * @tparam ValueType The type of entries to reduce.
	 */
	template< typename Operation, typename ValueType >
	class Reduce : public virtual MultiBSP_template {
	
		template< typename A, typename B >
		friend class AllReduce;

		private:

			/** Field that keeps track of the local reduction value. */
			ValueType reduced;

			/** Array of objects to reduce (locally, at leaf level only). */
			const ValueType *__restrict__ to_reduce;

			/** Size of the array to_reduce. */
			size_t length;

			/** Keeps track of which PID is root. */
			size_t target;

			/**
			 * Inner communication kernel of the all-to-one reduce,
			 * using bsp_direct_get.
			 *
			 * @param k Which sibling process to reduce from.
			 */
			void inner_all_to_one_reduce( const size_t k ) {
				//buffer used for direct_get
				ValueType buffer;
				//get remote contribution
				bsp_direct_get( k, &reduced, 0, &buffer, 1 );
				//reduce against local contribution
				Operation::apply( reduced, buffer );
			}

			/**
			 * Codes the communication kernel. Performs an all-to-one with
			 * interleaved reductions. Uses bsp_direct_get.
			 */
			void all_to_one_reduce() {
				//make sure we are ready to communicate
				bsp_sync();
				//get local PID
				const unsigned int s = bsp_pid();
				//if we are the root
				if( s == 0 ) {
					//do all-to-one + reduction here
					for( size_t k = 1; k < bsp_nprocs(); ++k ) {
						//reduce from PID k
						inner_all_to_one_reduce( k );
					}
					//check if we have target messages
					const void * payload;
					const void * tag;
					if( bsp_hpmove( &payload, &tag ) != SIZE_MAX ) {
						//we got a message! Now check for multiple roots
						if( bsp_hpmove( &payload, &tag ) != SIZE_MAX ) {
							bsp_abort( "Error in mcbsp::ops::Reduce, user defined more than one reduction root!\n" );
						}
						//set target
						target = true;
					}
				}
			}

			/** Retrieves global accumulation, if our MultiBSP subtree contains the root leaf process. */
			void get_accumulated() {
				//if we are root
				if( target ) {
					//get contribution from PID 0
					bsp_direct_get( 0, &reduced, 0, &reduced );
				}
			}

			/** Performs a local reduction of the array. */
			void localReduce() {
				//initialise result, use neutral element (i.e., the identity; which under
				//regular addition equals zero).
				reduced = Operation::identity;
				//loop over all input elements
				for( size_t i = 0; i < length; ++i ) {
					//reduce current element into result
					Operation::apply( reduced, to_reduce[ i ] );
				}
			}

		protected:

			/**
			 * Leaf-level code.
			 * If no leaf process is assigned root, the MultiBSP root program will be assigned
			 * root automatically. If more than one leaf process is set root, the program will
			 * abort.
			 *
			 * @param array The array to reduce.
			 * @param size  The length of the to-be reduced array.
			 * @param root  Whether this is the root process (false by default).
			 */
			ValueType leaf( const ValueType * const array, const size_t size, bool root = false ) {
				//set values
				to_reduce = array;
				length    = size;
				target    = root;
				//if we are root
				if( target ) {
					//let PID 0 know we are the target
					bsp_hpsend( 0, NULL, &target );
				}
				//do local reduce
				localReduce();
				//do communicate
				all_to_one_reduce();
				//let ancestors communicate, they will pass down the global result in PID 0
				bsp_up();
				//retrieve from PID 0, if necessary
				get_accumulated();
				//done
				return reduced;
			}

			/** Nested-level code. */
			ValueType nested() {
				//do all-to-one, reduces into PID 0
				all_to_one_reduce();
				//let ancestors communicate along PID 0
				bsp_up();
				//retrieve from PID 0, if necessary
				get_accumulated();
				//go back down to continue original execution flow
				bsp_down();
				//done
				return reduced;
			}

			/** Initialises the collective reduce operation. */
			void init() {
				bsp_push_reg( &reduced );
			}

			/** Releases the resources used for the collective reduce operation. */
			void destroy() {
				bsp_pop_reg( &reduced );
			}

			/**
			 * Reduce operation dispatcher. Automatically detects leaf or internal case. In
			 * case of a leaf call, the local array (or value) to be reduced is passed by a
			 * pointer to said array. Passing this array is optional; if none is passed, an
			 * empty array will be assumed. The globally reduced (accumulated) value will
			 * be returned by the root leaf MultiBSP program, as indicated by the root
			 * paramater to this call.
			 * If no leaf process is assigned root, the MultiBSP root program will be assigned
			 * as the reduction root automatically. If more than one leaf process is set root,
			 * the program will abort.
			 *
			 * Note that the array and size parameters have no meaning when calling the reduce
			 * operation from an internal MultiBSP program node; any non-NULL and non-zero
			 * values that may have been set will be ignored in the non-leaf case.
			 *
			 * @param array The array to reduce.
			 * @param size  The length of the to-be reduced array.
			 * @param root  Whether this is the root process (false by default).
			 *
			 * @return The global accumulated value, if root is set true. Otherwise, will
			 *         return a local accumulation only.
			 */
			ValueType apply( const ValueType * const array = NULL, const size_t size = 0, bool root = false ) {
				//leaf detection
				if( bsp_leaf() ) {
					//we are leaf node; call correct subprogram
					leaf( array, size, root );
				} else {
					//we are internal node; call correct subprogram
					nested();
				}
				//return global accumulation if we are root, otherwise return local accumulation
				return reduced;
			}

		public:

	};

} //collectives
} //mcbsp

#endif

