/*
 * Copyright (c) 2013 by Albert-Jan N. Yzelman
 *
 * This file is part of MulticoreBSP in C --
 *        a port of the original Java-based MulticoreBSP.
 *
 * MulticoreBSP for C is distributed as part of the original
 * MulticoreBSP and is free software: you can redistribute
 * it and/or modify it under the terms of the GNU Lesser 
 * General Public License as published by the Free Software 
 * Foundation, either version 3 of the License, or 
 * (at your option) any later version.
 * MulticoreBSP is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the 
 * implied warranty of MERCHANTABILITY or FITNESS FOR A 
 * PARTICULAR PURPOSE. See the GNU Lesser General Public 
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General 
 * Public License along with MulticoreBSP. If not, see
 * <http://www.gnu.org/licenses/>.
 */

/*! \file
 * Enables C++-style BSP communication; sizes and offsets
 * are given in number of elements of a specific type.
 * If any of the data pointers are of type (void*), then 
 * offset and sizes must be given in bytes as in plain C.
 */

#ifndef _HPP_MCBSP_TEMPLATES
#define _HPP_MCBSP_TEMPLATES

#ifndef _HPP_MCBSP
 #include "mcbsp.hpp"
#endif

/**
 * Registers a memory area for communication.
 *
 * This is the templated variant of a regular
 * bsp_push_reg().
 *
 * If an SPMD program defines a local variable x,
 * each of the P threads actually has its own memory
 * areas associated with that variable.
 * Communication requires threads to be aware of the
 * memory location of a destination variable. This
 * function achieves this.
 * The order of variable registration must be the
 * same across all threads in the SPMD program. The
 * size of the registered memory block may differ
 * from thread to thread. Registration takes effect
 * only after a synchronisation.
 *
 * Note: if MCBSP_COMPATIBILITY_MODE is defined,
 *       size will be of type `int'. Otherwise, it
 *       is of type `size_t'.
 *
 * @param address Pointer to the memory area to register.
 * @param size    Number of elements to register. This
 *                should equal 1 of address points to a
 *                single element instead of an array.
 */
template< typename T >
void bsp_push_reg( T * const address, const size_t size = 1 ) {
	bsp_push_reg( static_cast< void * >(address), size * sizeof( T ) );
}

/**
 * Sets the tag size of inter-thread messages.
 *
 * This is the templated variant of a regular
 * bsp_set_tagsize().
 *
 * bsp_send can be used to send message tuples
 * (tag,payload). Tag must be of a fixed size,
 * the payload size may differ per message.
 *
 * This function sets the tagsize so that it
 * can store an instance of a type equal to
 * the one passed as an argument.
 * This new size will be valid from the
 * next bsp_sync() on.
 *
 * All processes must still call 
 * bsp_set_tagsize with objects of the same
 * size, or MulticoreBSP for C will return an
 * error.
 *
 * The function will now output the old tag
 * size as the return value.
 *
 * Note: if MCBSP_COMPATIBILITY_MODE is defined,
 *       the return type is an `int'. Otherwise,
 *       it is of type `size_t'.
 *
 * @param tag_type An object of the new tag type.
 * @return The old tag size, in bytes.
 */

template< typename T >
MCBSP_BYTESIZE_TYPE bsp_set_tagsize( const T &tag_type = T() ) {
	MCBSP_BYTESIZE_TYPE size = static_cast< MCBSP_BYTESIZE_TYPE >( sizeof( tag_type ) );
	bsp_set_tagsize( &size );
	return size;
}

/**
 * Put data in a remote memory location.
 *
 * This is the templated variant of a regular bsp_put().
 *
 * This is a non-blocking communication request.
 * Communication will be executed during the next 
 * synchronisation step. The remote memory location must 
 * be registered using bsp_push_reg in a previous 
 * superstep.
 *
 * The data to be communicated to the remote area will 
 * be buffered on request; i.e., the source memory 
 * location is free to change after this communication 
 * request; the communicated data will not reflect 
 * those changes.
 *
 * Note: if MCBSP_COMPATIBILITY_MODE is defined, pid, 
 *       offset and size are of type `int'. Otherwise,
 *       pid is of type `unsigned int', and offset and
 *       size of type `size_t'.
 *
 * @param pid         The ID number of the remote thread.
 * @param source      Pointer to one or
 *                    more source elements.
 * @param destination Pointer to the registered remote 
 *                    memory area to send one or more
 *                    source elements to.
 * @param offset      Offset (in number of elements) of
 *                    the destination memory area.
 *                    Offset must be positive and must
 *                    not exceed the registered capacity
 *                    of the remote memory area.
 * @param size        Number of data elements to be
 *                    transmitted. I.e., all the data 
 *                    from source[0] up to source[size-1] 
 *                    at this process is copied into
 *                    destination[offset] up to
 *                    destination[offset+size-1] at the 
 *                    process with ID pid.
 */
template< typename T >
void bsp_put( const MCBSP_PROCESSOR_INDEX_DATATYPE pid, const void * const source,
	T * const destination, const size_t offset = 0, const size_t size = 1 ) {
	bsp_put( pid, source, static_cast< void * >(destination), offset * sizeof(T), size * sizeof(T) );
}

/**
 * Get data from a remote memory location.
 *
 * This is the templated variant of a regular bsp_get().
 *
 * This is a non-blocking communication request.
 * Communication will be executed during the next
 * synchronisation step. The remote memory location must
 * be registered using bsp_push_reg in a previous
 * superstep.
 *
 * The data retrieved will be the data at the remote
 * memory location at the time of synchronisation. It
 * will not (and cannot) retrieve data at "this"
 * point in the SPMD program at the remote thread.
 * If other communication at the remote process
 * would change the data at the region of interest,
 * these changes are not included in the retrieved
 * data; in this sense, the get is buffered.
 *
 * Note: if MCBSP_COMPATIBILITY_MODE is defined, pid, 
 *       offset and size are of type `int'. Otherwise,
 *       pid will be of type `unsigned int' and
 *       offset and size of type `size_t'.
 *
 * @param pid         The ID number of the remote thread.
 * @param source      Pointer to the registered remote 
 *                    memory area where to get data
 *                    elements from.
 * @param offset      Offset (in number of elements) of
 *                    the remote memory area.
 *                    Offset must be positive and must
 *                    not exceed the registered capacity
 *                    of the remote memory area.
 * @param destination Pointer to one or more local
 *                    destination elements.
 * @param size        Number of data elements to be
 *                    sent. I.e., all the data from
 *                    source[offset] up to
 *                    source[offset+size-1] at the
 *                    remote process is copied into 
 *                    destination[ 0 ] up to 
 *                    destination[size-1] at this
 *                    process.
 */
template< typename T >
void bsp_get( const MCBSP_PROCESSOR_INDEX_DATATYPE pid, const T * const source,
	const size_t offset, void * const destination, const size_t size = 1 ) {
	bsp_get( pid, static_cast< const void * >(source), offset * sizeof(T), destination, size * sizeof(T) );
}

/**
 * Get data from a remote memory location.
 *
 * This is the templated variant of a regular
 * bsp_direct_get().
 *
 * This is a blocking communication primitive:
 * communication is executed immediately and is not
 * queued until the next synchronisation step. The
 * remote memory location must be regustered using 
 * bsp_push_reg in a previous superstep.
 *
 * The data retrieved will be the data at the remote
 * memory location at "this" time. There is no
 * guarantee that the remote thread is at the same
 * position in executing the SPMD program; it might
 * be anywhere in the current superstep. If the 
 * remote thread writes to the source memory block
 * in this superstep, the retrieved data may
 * partially consist of old and new data; this
 * function does not buffer nor is it atomic in
 * any way.
 *
 * Note: if MCBSP_COMPATIBILITY_MODE is defined, pid, 
 *       offset and size are of type `int'. Otherwise,
 *       pid is of type `unsigned int' and offset &
 *       size of type `size_t'.
 *
 * @param pid         The ID number of the remote thread.
 * @param source      Pointer to the registered remote 
 *                    memory area where to get data
 *                    elements from.
 * @param offset      Offset (in number of elements) of
 *                    the remote memory area.
 *                    Offset must be positive and must
 *                    not exceed the registered capacity
 *                    of the remote memory area.
 * @param destination Pointer to one or more local
 *                    destination elements.
 * @param size        Number of data elements to be
 *                    sent. I.e., all the data from
 *                    source[offset] up to
 *                    source[offset+size-1] at the
 *                    remote process is copied into 
 *                    destination[ 0 ] up to 
 *                    destination[size-1] at this
 *                    process.
 */
template< typename T >
void bsp_direct_get( const MCBSP_PROCESSOR_INDEX_DATATYPE pid, const T * const source,
        const MCBSP_BYTESIZE_TYPE offset, void * const destination, const MCBSP_BYTESIZE_TYPE size = 1 ) {
	bsp_direct_get( pid, static_cast< const void * >(source), offset * sizeof(T), destination, size * sizeof(T) );
}

/**
 * Sends a message to a remote thread.
 *
 * This is the templated variant of a regular bsp_send().
 *
 * A message is actually a tuple (tag,payload). Tag is of a fixed size
 * (see bsp_set_tagsize), the payload size is set per message. Messages
 * will be available at the destination thread in the next superstep.
 *
 * Note: If MCBSP_COMPATIBILITY_MODE is defined, then pid and size are
 *       of type `int'. Otherwise, pid is of type `unsigned int' and
 *       size of type `size_t'.
 *
 * @param pid     ID of the remote thread to send this message to.
 * @param tag     Pointer to the tag data.
 * @param payload Pointer to one or more payload data elements.
 * @param size    Number of data elements in the payload.
 */
template< typename T >
void bsp_send( const MCBSP_PROCESSOR_INDEX_DATATYPE pid, const void * const tag,
	const T * const payload, const size_t size = 1 ) {
	bsp_send( pid, tag, static_cast< const void * >(payload), size * sizeof(T) );
}

#ifdef MCBSP_ENABLE_HP_DIRECTIVES
/**
 * This is a non-buffering and non-blocking send request.
 *
 * This is the templated variant of a regular bsp_hpsend().
 *
 * The function differs from the regular bsp_send in two major ways:
 * (1) the actual send communication may occur between now and the
 *     end of the next synchronisation phase, and
 * (2) the tag and payload data is read and sent somewhere between
 *     now and the end of the next synchronisation phase.
 * If you change the contents of the memory area tag and payload
 * point to after calling this function, undefined communication
 * will occur.
 * The semantics of BSMP remain unchanged: the sent messages will
 * only become available at the remote processor when the next
 * computation superstep begins. The performance gain is two-fold:
 * (a) bsp_hpsend avoids buffering-on-send, and
 * (b) BSP may send messages during a computation phase, thus
 *     overlapping computation with communication.
 *
 * Ad (a): normally, BSMP in BSP copies tag and payload data at
 * least three times. It buffers on bsp_send (buffer-on-send),
 * it buffers at the receiving processes' incoming BSMP queue 
 * (buffer-on-receive), and finally bsp_get_tag and bsp_move 
 * copy the data in the target user-supplied memory areas. To 
 * also eliminate the latter data movement, please consider
 * using bsp_hpmove.
 *
 * Ad (b): it is not guaranteed this overlap results in faster 
 * execution time. You should think about if using these high-
 * performance primitives makes sense on a per-application 
 * basis, and factor in the extra costs  of structuring your 
 * algorithm to enable correct use of these primitives.
 *
 * See bsp_send for general remarks about using BSMP primitives
 * in BSP. See bsp_hpget and bsp_hpput for equivalent (non-
 * buffering and non-blocking) high-performance communication
 * requests.
 *
 * Note: If MCBSP_COMPATIBILITY_MODE is defined, then pid and
 *       size are of type `int'. Otherwise, pid is of type 
 *       `unsigned int' and size of type `size_t'.
 *
 * @param pid     ID of the remote thread to send this message to.
 * @param tag     Pointer to the tag data. This data will not be
 *                buffered.
 * @param payload Pointer to one or more payload data elements.
 *                These elements will not be buffered.
 * @param size    Number of data elements in the payload.
 */
template< typename T >
void bsp_hpsend( const MCBSP_PROCESSOR_INDEX_DATATYPE pid, const void * const tag,
	const T * const payload, const size_t size = 1 ) {
	bsp_hpsend( pid, tag, static_cast< const void * >(payload), size * sizeof(T) );
}
#endif //MCBSP_ENABLE_HP_DIRECTIVES

/**
 * Retrieves the payload from the first message in the
 * queue of incoming messages, and removes that message.
 *
 * This is the templated variant of a regular bsp_move().
 *
 * If the incoming queue is empty, this function has
 * no effect. This function will copy a given maximum
 * of bytes from the message payload into a supplied
 * buffer. This maximum should equal or be larger
 * than the payload size (which can, e.g.,  be 
 * retrieved via bsp_get_tag).
 * The maximum can be 0 bytes; the net effect is
 * the efficient removal of the first message from
 * the queue.
 *
 * Note that Bulk Synchronous Message Passing (BSMP)
 * is doubly buffered: bsp_send buffers on send and
 * this function buffers again on receives.
 *
 * See bsp_hpmove if buffer-on-receive is unwanted.
 *
 * Note: if MCBSP_COMPATIBILITY_MODE is defined,
 *       max_copy_size is of type `int'. Otherwise,
 *       it is of type `size_t'.
 *
 * @param payload       Where to copy the data
 *                      elements from the payload of
 *                      a received BSMP message into.
 * @param max_copy_size The maximum number of
 *                      elements to copy.
 */
template< typename T >
void bsp_move( T * const payload, const size_t max_copy_size = 1 ) {
	bsp_move( static_cast< void * >(payload), max_copy_size * sizeof(T) );
}

//NOTE: there is no templated bsp_hpmove, as its intended use
//      requires casting the retrieved pointers into the BMSP
//      receive-buffer into the appropriate T*.
//      Additionally, the type of T is typically only known
//      after inspecting the tag; which is retrieved during
//      the same call to bsp_hpmove.

#if defined ENABLE_FAKE_HP_DIRECTIVES || defined MCBSP_ENABLE_HP_DIRECTIVES
/**
 * Put data in a remote memory location.
 *
 * This is the templated variant of a regular
 * bsp_hpput().
 *
 * Note: current implementation does a normal bsp_get!
 *       An communication overlapping implementation
 *       is forthcoming.
 *
 * This is a non-blocking communication request.
 * Communication will be executed sometime between now 
 * and during the next synchronisation step. Note that 
 * this differs from bsp_put. Communication is guaranteed 
 * to have finished before the next superstep. Note this
 * means that both source and destination memory areas 
 * might be read and written to at any time after 
 * issueing this communication request. This overlap of 
 * communication and computation is the fundamental
 * difference with the standard bsp_put.
 *
 * It is not guaranteed this overlap results in faster 
 * execution time. You should think about if using these 
 * high-performance primitives makes sense on a 
 * per-application basis, and factor in the extra costs 
 * of structuring your algorithm to enable correct use 
 * of these primitives.
 *
 * Otherwise usage is similar to that of bsp_put;
 * please refer to that function for further
 * documentation.
 *
 * Note: if MCBSP_COMPATIBILITY_MODE is defined, pid, 
 *       offset and size are of type `int'. Otherwise,
 *       pid is of type `unsigned int', offset and
 *       size of type `size_t'.
 *
 * @param pid         The ID number of the remote thread.
 * @param source      Pointer to one or
 *                    more source elements.
 * @param destination Pointer to the registered remote 
 *                    memory area to send one or more
 *                    source elements to.
 * @param offset      Offset (in number of elements) of
 *                    the destination memory area.
 *                    Offset must be positive and must
 *                    not exceed the registered capacity
 *                    of the remote memory area.
 * @param size        Number of data elements to be
 *                    transmitted. I.e., all the data 
 *                    from source[0] up to source[size-1] 
 *                    at this process is copied into
 *                    destination[offset] up to
 *                    destination[offset+size-1] at the 
 *                    process with ID pid.
 */
template< typename T >
void bsp_hpput( const MCBSP_PROCESSOR_INDEX_DATATYPE pid, const void * const source,
        T * const destination, const size_t offset = 0, const size_t size = 1 ) {
	bsp_hpput( pid, source, static_cast< void * >(destination), offset * sizeof(T), size * sizeof(T) );
}

/**
 * Get data from a remote memory location.
 *
 * This is the templated variant of a regular
 * bsp_hpget().
 *
 * Note: current implementation does a normal bsp_get!
 *       An communication overlapping implementation
 *       is forthcoming.
 *
 * This is a non-blocking communication request.
 * Communication will be executed between now and the
 * next synchronisation step. Note that this differs
 * from bsp_get. Communication is guaranteed to have
 * finished before the next superstep. Note this 
 * means that both source and destination memory 
 * areas might be read and written to at any time
 * after issueing this communication request. This
 * overlap of communication and computation is the
 * fundamental difference with the standard bsp_get.
 *
 * It is not guaranteed this overlap results in faster 
 * execution time. You should think about if using these 
 * high-performance primitives makes sense on a 
 * per-application basis, and factor in the extra costs 
 * of structuring your algorithm to enable correct use 
 * of these primitives.
 *
 * Note the difference between this high-performance
 * get and bsp_direct_get is that the latter function
 * is blocking (performs the communication
 * immediately and waits for it to end).
 *
 * Otherwise usage is similar to that of bsp_get;
 * please refer to that function for further
 * documentation.
 *
 * Note: if MCBSP_COMPATIBILITY_MODE is defined, pid, 
 *       offset and size are of type `int'. Otherwise,
 *       pid is of type `unsigned int' and offset &
 *       size of type `size_t'.
 *
 * @param pid         The ID number of the remote thread.
 * @param source      Pointer to the registered remote 
 *                    memory area where to get data
 *                    elements from.
 * @param offset      Offset (in number of elements) of
 *                    the remote memory area.
 *                    Offset must be positive and must
 *                    not exceed the registered capacity
 *                    of the remote memory area.
 * @param destination Pointer to one or more local
 *                    destination elements.
 * @param size        Number of data elements to be
 *                    sent. I.e., all the data from
 *                    source[offset] up to
 *                    source[offset+size-1] at the
 *                    remote process is copied into 
 *                    destination[ 0 ] up to 
 *                    destination[size-1] at this
 *                    process.
 */
template< typename T >
void bsp_hpget( const MCBSP_PROCESSOR_INDEX_DATATYPE pid, const T * const source,
        const size_t offset, void * const destination, const size_t size ) {
	bsp_hpget( pid, static_cast< void * >(source), offset * sizeof(T), destination, size * sizeof(T) );
}

#endif //ENABLE_FAKE_HP_DIRECTIVES

#endif

