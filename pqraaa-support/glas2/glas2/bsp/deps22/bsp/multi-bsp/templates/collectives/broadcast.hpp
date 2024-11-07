
#include "mcbsp-templates.hpp"
#include "multibsp.hpp"

#ifndef _H_MCBSP_BCAST
#define _H_MCBSP_BCAST

namespace mcbsp {
namespace collectives {

	/**
	 * Codes a broadcast collective operation. Here, a single leaf process, called the root,
	 * will send one data element to all other MultiBSP programs.
	 *
	 * @tparam ValueType The type of element to broadcast.
	 */
	template< typename ValueType >
	class Broadcast : public virtual MultiBSP_template {

		template< typename A, typename B >
		friend class AllReduce;

		private:

			/** Where the data from the root process will arrive at. */
			ValueType * receive_buffer;

			/** The buffer capacity. */
			size_t capacity;

			/** The number of ValueType elements to send within one broadcast operation. */
			size_t size;

			/** Whether we are the root of the communication. */
			bool root;

			/** Sends data to PID 0. */
			void send0() {
				bsp_send( 0, NULL, receive_buffer, size );
			}

			/** Executes neighbour communication, PID 0 to all others. */
			void bcast( const ValueType * const &toBroadcast ) {
				const unsigned int s = bsp_pid();
				const unsigned int P = bsp_nprocs();
				if( s == 0 ) {
					for( size_t k = 1; k < P; ++k ) {
						bsp_put( k, toBroadcast, receive_buffer, 0, size );
					}
				}
			}

			/** Sends root contribution to PID 0, and further upwards. */
			void up() {
				const unsigned int s = bsp_pid();
				//if I am root
				if( root ) {
					//if I am not PID 0
					if( s > 0 ) {
						//send contribution to PID 0
						send0();
					}
				}
				//allow PID 0 to receive broadcast value
				bsp_sync();
				//if I am PID 0
				if( s == 0 ) {
					//check if I received messages
					unsigned int num_msg;
					bsp_qsize( &num_msg, NULL );
					if( num_msg > 0 ) {
						//if I am root I should not have received anything
						if( root ) {
							bsp_abort( "Error: multiple root processes in Broadcast operation!\n" );
						} else {
							//I should never receive more than 1 contribution
							if( num_msg > 1 ) {
								bsp_abort( "Error: multiple root processes in Broadcast operation!\n" );
							}
						}
						//all OK, received the broadcast element
						bsp_move( receive_buffer, size );
						root = true;
					}
				}
				//allow broadcast element to bubble upwards
				bsp_up();
			}

		public:

			/** Default constructor. */
			Broadcast() : receive_buffer( NULL ), capacity( 0 ) {}

			/** Initialises this instance for collective broadcast operations. */
			void init( const size_t _cap = 1 ) {
				//initialisation here is slightly involved. We have to make sure
				//to execute the leaf-level initialisation first. TODO
				//sanity check
				assert( _cap > 0 );
				//check if action is required at all
				if( capacity >= _cap ) {
					//current buffer capacity suffices, register it
					bsp_push_reg( receive_buffer, capacity );
					//and exit
					return;
				}
				//new buffer must be allocated; first destroy existing buffer, if there is one
				if( capacity > 0 ) {
					destroy();
				}
				//leaf processes allocate buffer
				if( bsp_leaf() ) {
					//set new capacity
					capacity = _cap;
					//allocate new buffer
					receive_buffer = new ValueType[ capacity ];
				}
				//sanity check
				assert( receive_buffer != NULL );
				//DBG
				//std::cout << "Registered receive buffer @ " << receive_buffer << " leaf status = " << bsp_leaf() << std::endl;
				//register new buffer
				bsp_push_reg( receive_buffer, capacity );
				//done
			}

			/** Undo all initialisation performed for collective broadcast operations. */
			void destroy() {
				//pop registry
				//bsp_pop_reg( receive_buffer );
				//leaf processes destroy buffer
				if( bsp_leaf() ) {
					delete [] receive_buffer;
				}
			}

			/**
			 * Broadcast code on the leaf function.
			 * @param toBroadcast When not NULL, this is a pointer to the value to be
			 *        broadcasted, and the calling process is the root of this broadcast 
			 *        operation. When NULL, this process will participate in this 
			 *        collective operation as a non-root process.
			 * @param _s The length of the to-be broadcasted vector.
			 * @return The broadcasted value.
			 */
			const ValueType * leaf( const size_t _s ) {
				return leaf( NULL, _s );
			}

			/**
			 * Broadcast code on the leaf function.
			 * @param toBroadcast When not NULL, this is a pointer to the value to be
			 *        broadcasted, and the calling process is the root of this broadcast 
			 *        operation. When NULL, this process will participate in this 
			 *        collective operation as a non-root process.
			 * @param _s The length of the to-be broadcasted vector.
			 * @return The broadcasted value.
			 */
			const ValueType * leaf( const ValueType * const toBroadcast, const size_t _s ) {
				//sanity check
				if( _s > capacity ) {
					bsp_abort( "To-be broadcasted vector is of larger length than the buffer capacity! Please call Broadcast::init with a more appropriate argument!\n" );
				}
				//set size of this operation
				size = _s;
				//check for root
				if( toBroadcast == NULL ) {
					//no root
					root = false;
				} else {
					//set to-be broadcasted values
					for( size_t i = 0; i < size; ++i ) {
						receive_buffer[ i ] = toBroadcast[ i ];
					}
					//we are root
					root = true;
				}
				//root contribution needs to be sent to the MultiBSP root
				up();
				//PID 0 sends contribution to siblings
				bcast( receive_buffer );
				//allow communication to finish
				bsp_sync();
				//done; return broadcasted value
				return receive_buffer;
			}

			/**
			 * Broadcast code on the nested level.
			 * @param child Reference to the nested PID 0 data (obtainable via bsp_retrieve).
			 * @param _s The length of the to-be broadcasted vector.
			 * @return The broadcasted value.
			 */
			void nested() {
				//sanity check
				if( size > capacity ) {
					bsp_abort( "A leaf process has set Broadcast::size to a larger length than the buffer capacity; please call Broadcast::init with a more appropriate argument!\n" );
				}
				//propagate upwards
				up();
				//upstream has set the broadcasted value, distribute it amongst our siblings
				bcast( receive_buffer );
				//finish comm and go down
				bsp_down();
				//done
			}

	};

} //collectives
} //mcbsp

#endif

