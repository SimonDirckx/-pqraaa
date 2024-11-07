/*
 * Copyright (c) 2007-2014, A. N. Yzelman,   Utrecht University 2007-2011;
 *                                                    KU Leuven 2011-2014.
 *                          R. H. Bisseling, Utrecht University 2007-2014.
 * 
 * This file is part of the Sparse Library.
 * 
 * This library was developed under supervision of Prof. dr. Rob H. Bisseling at
 * Utrecht University, from 2007 until 2011. From 2011-2014, development continued 
 * at KU Leuven, where Prof. dr. Dirk Roose contributed significantly to the ideas 
 * behind the newer parts of the library code.
 * 
 *     The Sparse Library is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by the
 *     Free Software Foundation, either version 3 of the License, or (at your
 *     option) any later version.
 * 
 *     The Sparse Library is distributed in the hope that it will be useful, but
 *     WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 *     or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 *     for more details.
 * 
 *     You should have received a copy of the GNU General Public License along
 *     with the Sparse Library. If not, see <http://www.gnu.org/licenses/>.
 */

/* \file
 * File created by:
 *     A. N. Yzelman, Dept. of Computer Science, KU Leuven, 2011.
 */

#include <map>
#include <vector>
#include <iostream>
#include <limits.h>
#include <pthread.h>

#ifndef _NO_LIBNUMA
 #include <numa.h>
#endif

#ifndef SIZE_MAX
 #define SIZE_MAX static_cast< size_t >( -1 )
#endif

#include "SparseMatrix.hpp"
#include "Matrix2HilbertCoordinates.hpp"
#include "FBICRS.hpp"
#include "MachineInfo.hpp"

#ifndef _H_BETAHILBERT
#define _H_BETAHILBERT

/** 
 * Collection method for gathering local results and gathering them in the global y.
 *
 * 0: sequential get; 	one thread copies all local data to the global output vector
 * 1: parallel get;   	global output vector is block distributed and all threads do a get on their local part
 * 2: reduction-based; 	output vectors are block distributed, collection uses a reduction tree
 * 			(reduction base is defined in BH_REDUCE_BASE)
 * 3: unscalable red.; 	reduction which does not scale; no block distribution, distributing leafs
 *
 *-1: NO COLLECT: 	for debug purposes only!!!
 */
#define BH_COLLECT 1

/**
 *  Reduction base for collect method 2 (see above). On a dual dual-core machine,
 *  the proper base is two; on a dual quadcore machine, the proper base is either
 *  two or four; on a dual 10-core machine the proper choice is 10, and so on.
 */
#define BH_REDUCE_BASE 2

/**
 * Thread 0 will use BetaHilbert< T >::output as `local y'. This y is not allocated 
 * by that thread, however, and thus may incur a performance penalty if that output
 * vector was allocated on another socket, for instance. When you suspect this
 * happens, turn this paramater OFF.
 *
 * 0: off
 * 1: on
 */
#define BH_USE_GLOBAL_Y 1

/**
 *  Shared data for BetaHilbert threads.
 */
template< typename T >
class shared_data {

	public:

		/** Thread ID */
		size_t id;

		/** Total number of processors */
		size_t P;

		/** 0 undef, 1 init, 2 zax, 3 zxa, 4 exit, 5 reset */
		unsigned char mode;

		/** how many times to repeat the operation set in `mode' (above, only for 2 and 3) */
		size_t repeat;

		/** Array of local sparse blocks. */
		std::vector< Triplet< T > > *original;

		/** Will cache block numbers of nonzeroes */
		size_t *nzb;

		/** Will contain the nonzero counts of separate blocks */
		size_t **nzc;

		/** Will store local timing */
		double time;

		/** Local memory use */
		size_t bytes;

		/** Mutex used for synchronisation. */
		pthread_mutex_t* mutex;

		/** Condition used for synchronisation. */
		pthread_cond_t*  cond;

		/** Mutex used for end sync. */
		pthread_mutex_t* end_mutex;

		/** Condition used for end sync. */
		pthread_cond_t*  end_cond;

		/** Counter used for synchronisation. */
		size_t *sync;

		/** Counter used for end sync. */
		size_t *end_sync;

		/** Length of the local output vector. */
		size_t output_vector_size;

		/** Offset of the local output vector compared to global indices. */
		size_t output_vector_offset;

		/** Pointer to the local output vector. */
		T *local_y;

		/** Array of all local input vectors of all SPMD processes. */
		const T **input;

		/** Array of all output vectors local to all SPMD processes. */
		T **output;

		/** Base constructor. Initialises with invalid default arguments. */
		shared_data(): id( -1 ), P( -1 ), mode( 0 ), original( NULL ), nzb( NULL ), nzc( NULL ), time( 0 ),
				mutex( NULL ), cond( NULL ), end_mutex( NULL ), end_cond( NULL ), sync( NULL ),
				end_sync( NULL ), output_vector_size( -1 ), output_vector_offset( -1 ),
				local_y( NULL ), input( NULL ), output( NULL ) {}

		/**
		 * Default constructor.
		 *
		 * @param _id Thread ID.
		 * @param _P  Number of SPMD processes.
		 * @param _original Original set of nonzeroes.
		 * @param _nzb Number of sparse blocks.
		 * @param _nzc Number of nonzeroes in each sparse block,
		 * @param _mutex Sync mutex.
		 * @param _cond  Sync condition.
		 * @param _end_mutex End sync mutex.
		 * @param _end_cond   End sync condition.
		 * @param _sync Sync counter.
		 * @param _end_sync End sync counter.
		 * @param _m Number of matrix rows.
		 * @param _in Input vectors of all SPMD processes.
		 * @param _out Output vectors of all SPMD processes.
		 */
		shared_data( size_t _id, size_t _P, std::vector< Triplet< double > > *_original, size_t *_nzb, size_t **_nzc, 
				pthread_mutex_t *_mutex, pthread_cond_t *_cond, pthread_mutex_t *_end_mutex, pthread_cond_t *_end_cond,
				size_t *_sync, size_t *_end_sync, size_t _m,
				const T **_in, T**_out ):
				id( _id ),  P( _P ), mode( 1 ), original( _original ), nzb( _nzb ), nzc( _nzc ), time( 0 ),
				mutex( _mutex ), cond( _cond ), end_mutex( _end_mutex ), end_cond( _end_cond ),
				sync( _sync ), end_sync( _end_sync ), output_vector_size( _m ), output_vector_offset( -1 ),
				local_y( NULL ), input( _in ), output( _out ) {}
};

/** Compare function used for quicksort on an array of size_ts */
int beta_uli_compare( const void *a, const void *b ) {
	return ( *(size_t*)a - *(size_t*)b );
}

/** The Beta Hilbert triplet scheme. Full parallel SpMV, based on Blocked Hilbert and PThreads.
 *  Inspired by Aydin & Gilbert's CSB, and comments by Patrick Amestoy on the BICRS Hilbert scheme. */
template< typename T >
class BetaHilbert: public SparseMatrix< T, ULI > {

	private:

	protected:

		/** Number of threads to fire up */
		static size_t P;

		/** Which processors to pin threads to */
		std::vector< unsigned short int > p_translate;

		/** Input vector */
		const T* input;

		/** Output vector */
		T* output;

		/** p_threads associated to this data strcuture */
		pthread_t *threads;

		/** array of initial thread data */
		shared_data< T > *thread_data;

		/** Clock type used for thread-local timing */
		static clockid_t global_clock_id;

		/** Stop/continue mechanism: mutex */
		pthread_mutex_t mutex;

		/** Stop/continue mechanism: condition */
		pthread_cond_t cond;

		/** Wait for end mechanism: mutex */
		pthread_mutex_t end_mutex;

		/** Wait for end mechanism: condition */
		pthread_cond_t end_cond;

		/** Used for synchronising threads */
		size_t sync;

		/** Used for construction end signal */
		size_t end_sync;

		/** Given FBICRS, the maximum value for columnwise matrix size */
		static const ULI max_n = FBICRS< T >::beta_n;

		/** Given FBICRS, the maximum value for the rowwise matrix size, assuming short ints on ICRS at the lower level */
		static const ULI max_m = FBICRS< T >::beta_m;

		/** Sets p_translate to 0..P-1 by default, or equal to the optionally supplied vector. */
		void set_p_translate( std::vector< unsigned short int > *_p_translate ) {
			if( _p_translate == NULL ) {
				//get number of cores available
				P = MachineInfo::getInstance().cores();
				for( unsigned short int i = 0; i < P; ++i ) p_translate.push_back( i );
			} else {
				p_translate = *_p_translate;
				P = p_translate.size();
			}
		}

	public:

		/**
		 * File-based constructor.
		 * 
		 * @param file Filename to read from.
		 * @param zero Zero element of the sparse matrix.
		 * @param _p_translate Thread pinning array.
		 */
		BetaHilbert( const std::string file, T zero = 0, std::vector< unsigned short int > *_p_translate = NULL ): input( NULL ), output( NULL ) {
			set_p_translate( _p_translate );
			this->loadFromFile( file, zero );
		}

		/**
		 * COO-based constructor.
		 * 
		 * @param input COO matrix.
		 * @param m The number of rows in the input matrix.
		 * @param n The number of columns in the input matrix.
		 * @param zero Zero element of the sparse matrix.
		 * @param _p_translate Thread pinning array.
		 */
		BetaHilbert( std::vector< Triplet< T > >& input, ULI m, ULI n, T zero = 0, std::vector< unsigned short int > *_p_translate = NULL ): input( NULL ), output( NULL ) {
			set_p_translate( _p_translate );
			load( input, m, n, zero );
		}

		/** Base deconstructor. */
		virtual ~BetaHilbert() {
			//set all daemon threads to exit mode
			for( size_t i=0; i<P; i++ ) {
				thread_data[ i ].mode = 4;
			}

			//wake up all daemon threads
			pthread_mutex_lock( &mutex );
			pthread_cond_broadcast( &cond );
			pthread_mutex_unlock( &mutex );

			//allow threads to exit gracefully
			for( size_t i=0; i<P; i++ )
				pthread_join( threads[ i ], NULL );

			//destroy data
			delete [] thread_data;
			delete [] threads;
			pthread_mutex_destroy( &mutex );
			pthread_cond_destroy(  &cond  );
		}

		/** Callee will Wait for end of SpMV. */
		void wait() {
			//wait for end signal
			pthread_cond_wait( &end_cond, &end_mutex );
			pthread_mutex_unlock( &end_mutex );
		}

		/**
		 * Loads from input COO matrix.
		 *
		 * @param input The input COO matrix.
		 * @param m Input matrix row dimension.
		 * @param n Input matrix column dimension.
		 * @param zero Zero element of the input matrix.
		 */
		virtual void load( std::vector< Triplet< T > >& input, const ULI m, const ULI n, const T zero ) {
#ifndef _NO_LIBNUMA
			//set kernel to local thread allocation if it wasn't already the case
			numa_set_localalloc();
#endif

			//base settings
			this->zero_element = zero;
			this->nor = m;
			this->noc = n;
			this->nnz = input.size();

			// Parallel distribution according to the Hilbert curve:
			// -get maximum Hilbert coordinate
			//  (Note: curve coordinates should be implented such that (0,0) maps to 0!
			//         Otherwise a minimum Hilbert coordinate should also be found, and
			//         below code must be adapted to work on a Hilbert coordinate range
			//         instead)
			// -simultaneously also get number of nonzeroes to locally store
			// -translate to a range
			// -call function to distribute this range over P threads
			//
			// Most operations are done in parallel (see thread function), basically
			// only shared array allocation and the actual forking is in this function.

			size_t  *nzb = new size_t [ this->nnz + P ];
			size_t **nzc = new size_t*[ P ];
			for( size_t i=0; i<P; i++ ) nzc[ i ] = NULL;

			//create P threads :)
			this->threads = new pthread_t[ P ];
			//initialize local initialisation data
			thread_data = new shared_data< T >[ P ];
			//initialize mutexes and conditions and a synchronisation counter
			pthread_mutex_init( &mutex, NULL );
			pthread_cond_init ( &cond,  NULL );
			pthread_mutex_init( &end_mutex, NULL );
			pthread_cond_init ( &end_cond,  NULL );
			sync     = 0;
			end_sync = 0;
			//lock end mutex (disallow threads that are created to signal for end
			//before this thread is done with spawning children)
			pthread_mutex_lock( &end_mutex );
			//go forth and multiply
			for( size_t i=0; i<P; i++ ) {
				//build thread-local init data
				thread_data[ i ] = shared_data< T >( i, P, &input, nzb, nzc, &mutex, &cond, &end_mutex, &end_cond, &sync, &end_sync, this->m(), &(this->input), &output );
				//set fixed affinity for threads
				cpu_set_t mask;
				CPU_ZERO( &mask );
				CPU_SET ( p_translate[ i ], &mask );
				//prepare attributes
				pthread_attr_t attr;
				pthread_attr_init( &attr );
				//set fixed affinity in attribute, so that it starts binded immediately
				pthread_attr_setaffinity_np( &attr, sizeof( cpu_set_t ), &mask );
				//fire up thread
				pthread_create( &threads[ i ], &attr, &BetaHilbert::thread, (void*) &thread_data[i] );
				//free attr
				pthread_attr_destroy( &attr );
			}

			//wait for threads to finish initialisation
			wait();

			//delete temporary array
			delete [] nzb;
			delete [] nzc;
		}

		/**
		 * End synchronisation function.
		 *
		 * @param mutex Which mutex to sync with.
		 * @param cond  Which condition to sync with.
		 * @param sync  Sync counter to be used.
		 * @param P     How many SPMD processes are running.
		 */
		static void end( pthread_mutex_t* mutex, pthread_cond_t* cond, size_t *sync, const size_t P ) {
			pthread_mutex_lock( mutex );
			(*sync)++;
			if( *sync == P ) {
				pthread_cond_signal( cond );
				*sync = 0;
			}
			pthread_mutex_unlock( mutex );
		}

		/**
		 * Thread synchronisation function.
		 * 
		 * @param mutex Mutex used for synchronisation.
		 * @param cond  Condition used for synchronisation.
		 * @param sync  Synchronisation counter.
		 * @param P     The number of SPMD processes running.
		 */
		static void synchronise( pthread_mutex_t* mutex, pthread_cond_t* cond, size_t *sync, const size_t P ) {
			pthread_mutex_lock( mutex );
			(*sync)++;
			if( *sync == P ) {
				*sync = 0;
				pthread_cond_broadcast( cond );
			} else
				pthread_cond_wait( cond, mutex );
			pthread_mutex_unlock( mutex );
		}

		/**
		 * SPMD code for each thread doing an SpMV.
		 * 
		 * @param data Pointer to the thread-local data.
		 * @return Pointer to the thread output data.
		 */
		static void* thread( void *data ) {
			shared_data<T>* shared = (shared_data<T>*)data;
			const size_t id  = shared->id;
			const size_t P   = shared->P;
			const size_t nnz = shared->original->size();
			pthread_mutex_t *mutex = shared->mutex;
			pthread_cond_t  *cond  = shared->cond;

			//check whether pinned to exactly one unit of execution
			cpu_set_t mask;
			CPU_ZERO( &mask );
			pthread_getaffinity_np( pthread_self(), sizeof( cpu_set_t ), &mask );
			for( size_t s=0; s<P; s++ ) {
				if( s==id ) continue;
				if( CPU_ISSET( s, &mask ) ) {
					std::cerr << "Thread " << id << " mask is larger than one core" << " (" << s << " is set)!" << std::endl;
					exit( 1 );
				}
			}
#ifndef NDEBUG
			std::cout << "Phase 1 at thread " << shared->id << ": cache Hilbert values and get maximum one" << std::endl;
#endif
			size_t h1, h2, max1, max2;
			size_t *h2s = shared->nzb;
			const size_t blocksize = (nnz % P) > 0 ? nnz / P + 1 : nnz / P;
			Matrix2HilbertCoordinates::IntegerToHilbert( (*(shared->original))[ id*blocksize ].i() / max_m, (*(shared->original))[ id*blocksize ].j() / max_n, h1, h2 );
			h2s[ id * blocksize ] = h2;
			max1 = h1;
			max2 = h2;
			for( size_t i=id*blocksize+1; i<shared->original->size() && i < (id+1)*blocksize; ++i ) {
				Matrix2HilbertCoordinates::IntegerToHilbert( (*(shared->original))[ i ].i() / max_m,
					(*(shared->original))[ i ].j() / max_n, h1, h2 );
				h2s[ i ] = h2;
				if( h1 >  max1 ) { max1 = h1; max2 = h2; }
				else if( h1 == max1 && h2 > max2 ) max2 = h2;
			}

			//check for too large ranges
			if( max1 > 0 ) {
				std::cerr << "Hilbert coordinate range is larger than 2^64-1. Current counting sort mechanism thus requires an array of 2^64 integers-->quitting." << std::endl;
				//it is not possible to allocate such arrays on 64-bit machines. Also, it is extremely unlikely counting sort is the most efficient sort when this happens--
				//very likely most Beta blocks will be empty and a normal quicksort is more appropriate.
				exit( 1 );
			}

			//communicate local maxima
			h2s[ nnz + id ] = max2; //false sharing, but well...

			//sync
			BetaHilbert::synchronise( mutex, cond, shared->sync, shared->P );

			//determine global maximum
			for( size_t i=0; i<P; ++i )
				if( h2s[ shared->original->size() + i ] > max2 )
					max2 = h2s[ shared->original->size() + i ];
			max2++; //number of betablocks is maximum index plus 1

			//check for block size appropriateness
			unsigned char proceed;
			if( max2 >= nnz ) {
#ifndef NDEBUG
				std::cout << "Number of beta_m x beta_n blocks exceeds actual number of nonzeroes; going for (sequential) quicksort implementation." << std::endl;
#endif
				proceed = 0;
			} else {
#ifndef NDEBUG
				std::cout << "Choosing counting-sort implementation." << std::endl;
#endif
				proceed = 1;
			}

			// These values are to be filled
			size_t m, n, start, end;
			std::vector< std::vector< Triplet< T > > > beta;

			if( proceed == 0 ) {
				//temporary array
				size_t *horig = NULL;

				//sequential quicksort O(nz*log(nz))
				if( id == 0 ) {
					horig = new size_t[ shared->original->size() ];
					for( size_t i=0; i<shared->original->size(); ++i )
						horig[ i ] = h2s[ i ];
					qsort( h2s, nnz, sizeof( size_t ), beta_uli_compare );
					assert( h2s[ nnz-1 ]+1 == max2 );
				}
				//sync
				BetaHilbert::synchronise( mutex, cond, shared->sync, shared->P );
				//determine distribution
				/*for( size_t i = blocksize, k=0; i < shared->original->size(); i += blocksize, ++k ) {
					const size_t target = h2s[ i ];
					size_t split = i;
					for( size_t dev_f = 0, dev_b = 0; ; dev_f++, dev_b-- ) {
						if( i + dev_f < shared->original->size() ) {
							if( h2s[ i + dev_f ] != target ) {
								split = i + dev_f - 1;
								i += dev_f / 2;
								break;
							}
						} else if( i - dev_b >= 0  ) {
							if( h2s[ i - dev_b ] != target ) {
								split = i - dev_b;
								i += dev_b / 2;
								break;
							}
						} else
							break;
					}
					if( shared->id == k ) {
						std::cout << "Processor " << shared->id << " chooses to split at " << split << " with h=" << h2s[ split ] << " while target=" << target << std::endl;
						start = h2s[ split ];
					} else if( shared->id + 1 == k ) {
						std::cout << "Processor " << shared->id << " chooses to end at " << split << " with h=" << h2s[ split ] << "." << std::endl;
						end = h2s[ split ];
					}
				}
				if( shared->id == shared->P - 1 ) end = h2s[ shared->original->size() - 1 ];*/
				start = id * blocksize;
				end   = start + blocksize;
				if( end > max2 ) end = max2;
				const size_t target_s = h2s[ start ];
				while( start < max2 && h2s[ start ] == target_s ) start--;
				start++;
				const size_t target_e = h2s[ end ];
				while( end < max2 && h2s[ end ] == target_e ) end--;
				end++;
				if( shared->id + 1 == shared->P ) end = nnz - 1;
				start = h2s[ start ];
				end   = h2s[ end ];
#ifndef NDEBUG
				std::cout << "Processor " << id << " range is " << start << " to " << end << ", max = " << max2 << std::endl;
#endif
				assert( start <= end );
				assert( end <= max2 );
				
				//build local data structure from the triplets in shared->original
				//and beta block hilbert index to beta translation map
				std::map< size_t, size_t > h2b;

				//prepare beta and map
				//get first hilbert block coordinate
				h2 = h2s[ 0 ];
				//create corresponding vector in beta
				beta.push_back( std::vector< Triplet< T > >() );
				//remember index
				h2b[ h2 ] = 0;
				//keep next index
				size_t c = 1;
				//do the same for the remainder of available hilbert coordinates
				for( size_t i=1; i<shared->original->size(); ++i ) {
					h2 = h2s[ i ];
					beta.push_back( std::vector< Triplet< T > >() );
					h2b[ h2 ] = c++;
				}
					
				//prepare to get matrix size (m, n) and some statistics too
				size_t smin_m = SIZE_MAX;
				size_t smin_n = SIZE_MAX;
				size_t smax_m, smax_n;
				smax_m = smax_n = m = n = 0;
				size_t *ms    = new size_t[ end - start ];
				size_t *ns    = new size_t[ end - start ];
				size_t *minms = new size_t[ end - start ];
				size_t *minns = new size_t[ end - start ];
				for( size_t i=0; i<end-start; i++ ) {
					ms[ i ] = 0;
					ns[ i ] = 0;
					minms[ i ] = SIZE_MAX;
					minns[ i ] = SIZE_MAX;
				}

				if( id == 0 ) {
					for( size_t i=0; i<nnz; i++ )
						h2s[ i ] = horig[ i ];
					delete [] horig;
				}

				for( size_t i=0; i<shared->original->size(); i++ ) {
					const ULI row = (*(shared->original))[ i ].i();
					const ULI col = (*(shared->original))[ i ].j();
					h2 = h2s[ i ];
					if( row > m ) m = row;
					if( col > n ) n = col;
					if( h2 >= start && h2 < end ) {
						if( row > ms[ h2-start ] ) ms[ h2-start ] = row;
						if( col > ns[ h2-start ] ) ns[ h2-start ] = col;
						if( row < minms[ h2-start ] ) minms[ h2-start ] = row;
						if( col < minns[ h2-start ] ) minns[ h2-start ] = col;
						if( row < smin_m ) smin_m = row;
						if( col < smin_n ) smin_n = col;
						if( row > smax_m ) smax_m = row;
						if( col > smax_n ) smax_n = col;
					}
				}
				//matrix size is 1 plus max index
				m++;
				n++;

				//sync
				BetaHilbert::synchronise( mutex, cond, shared->sync, shared->P );

				for( size_t i=0; i<shared->original->size(); i++ ) {
					h2 = h2s[ i ];
					assert( h2 < max2 );
					if( h2 >= start && h2 < end ) {
						beta[ h2b[ h2 ] ].push_back( (*(shared->original))[ i ] );
					}
				}
				
				size_t cursum = 0;
				for( size_t i=0; i<beta.size(); i++ )
					cursum += beta[i].size();
#ifndef NDEBUG
				std::cout << "Thread " << shared->id << ": " << smin_m << "," << smax_m << " times " << smin_n << "," << smax_n << " holding " << cursum << " nonzeroes." << std::endl;
#endif

				//sync
				BetaHilbert::synchronise( mutex, cond, shared->sync, shared->P );

				//remove temporary values
				delete [] ms;
				delete [] ns;
				delete [] minms;
				delete [] minns;

			} else 	if( proceed == 1 ) {
				//start counting sort
#ifndef NDEBUG
				std::cout << "Phase 2: count nonzeroes in each beta_m x beta_n block (counting sort step 1), and derive distribution" << std::endl;
#endif
				size_t *nzc = new size_t[ max2 ];
				for( size_t i = 0; i < max2; ++i ) nzc[ i ] = 0;
				for( size_t i=id*blocksize; i<nnz && i < (id+1)*blocksize; ++i ) nzc[ h2s[ i ] ]++;
				shared->nzc[ id ] = nzc;
				
				//sync
				BetaHilbert::synchronise( mutex, cond, shared->sync, shared->P );
	
				size_t sum = 0;
#ifndef NDEBUG
				//sanity check
				if( id == 0 ) {
					for( size_t i=0; i<max2; i++ )
						for( size_t s=0; s<P; s++ )
							sum += shared->nzc[ s ][ i ];
					assert( sum == nnz );
					sum = 0;
				}
				BetaHilbert::synchronise( mutex, cond, shared->sync, shared->P );
#endif

				//combine counting sort
				for( size_t i=id*(max2/P+1); i<(id+1)*(max2/P+1)&&i<max2; ++i ) {
					for( size_t s=1; s<P; ++s )
						shared->nzc[ 0 ][ i ] += shared->nzc[ s ][ i ];
				}

				//sync
				BetaHilbert::synchronise( mutex, cond, shared->sync, shared->P );

#ifndef NDEBUG
				//sanity check
				sum = 0;
				for( size_t i=0; i<max2; i++ ) sum += shared->nzc[ 0 ][ i ];
				assert( sum == nnz );
				sum = 0;
#endif

				//get start and end positions of processor 0
				const size_t nnz_target = nnz != 0 ? nnz / P + 1 : nnz / P;
				size_t start  = 0;
				size_t end    = 1;
				size_t cursum = shared->nzc[ 0 ][ 0 ];
				for( ; end<max2 && cursum < nnz_target; end++ )
					cursum += shared->nzc[ 0 ][ end ];
			
				//get my start and end position
				for( size_t i=0; i<P; i++ ) {
					if( shared->id == i ) {
						std::cout << "Thread " << i << ": local nnz count is " << cursum << " (storing block indices " << start << ", " << end << ") out of " << max2 << " blocks present." << std::endl;
						break;
					}
					sum += cursum;
					start = end;
					cursum = 0;
					for( ; end < max2 && ( (i == P-1) || (cursum < nnz_target) ); end++ )
						cursum += shared->nzc[ 0 ][ end ];
				}

#ifndef NDEBUG
				std::cout << "Phase 3 at processor " << shared->id << ": getting local nonzeroes (counting sort step 2)" << std::endl;
#endif

				//build local data structure from the triplets in shared->original
				beta.resize( end - start );
				for( size_t i=0; i<end-start; i++ )
					beta[ i ] = std::vector< Triplet< T > >( shared->nzc[ 0 ][ i + start ] );

				//prepare to get matrix size (m, n) and some statistics too
				size_t smin_m, smin_n;
				size_t smax_m, smax_n;
				smin_m = smin_n = ULONG_MAX;
				smax_m = smax_n = m = n = 0;
				size_t *ms    = new size_t[ end - start ];
				size_t *ns    = new size_t[ end - start ];
				size_t *minms = new size_t[ end - start ];
				size_t *minns = new size_t[ end - start ];
				for( size_t i=0; i<end-start; i++ ) {
					ms[ i ] = 0;
					ns[ i ] = 0;
					minms[ i ] = ULONG_MAX;
					minns[ i ] = ULONG_MAX;
				}

				//copy local elements
				for( size_t i=0; i<end-start; i++ )
					assert( beta[ i ].size() == shared->nzc[ 0 ][ i + start ] );

				for( size_t i=0; i<nnz; ++i ) {
					const ULI row = (*(shared->original))[ i ].i();
					const ULI col = (*(shared->original))[ i ].j();
					h2 = shared->nzb[ i ];
					if( row > m ) m = row;
					if( col > n ) n = col;
					if( h2 >= start && h2 < end ) {
						if( row > ms[ h2-start ] ) ms[ h2-start ] = row;
						if( col > ns[ h2-start ] ) ns[ h2-start ] = col;
						if( row < minms[ h2-start ] ) minms[ h2-start ] = row;
						if( col < minns[ h2-start ] ) minns[ h2-start ] = col;
						if( row < smin_m ) smin_m = row;
						if( col < smin_n ) smin_n = col;
						if( row > smax_m ) smax_m = row;
						if( col > smax_n ) smax_n = col;
					}
				}
				//matrix size is 1 plus max index
				m++;
				n++;

				//sync here because the following phase will modify nzc,
				//while it may still be used above (and assumed constant)
				BetaHilbert::synchronise( mutex, cond, shared->sync, shared->P );

				for( size_t i=0; i<nnz; ++i ) {
#ifndef NDEBUG
					const ULI row = (*(shared->original))[ i ].i();
					const ULI col = (*(shared->original))[ i ].j();
#endif
					h2 = shared->nzb[ i ];
					assert( h2 < max2 );
					if( h2 >= start && h2 < end ) {
#ifndef NDEBUG
						//sanity check: this nonzero should be grouped with other nonzeroes in the same block
						if( beta[ h2-start ].size() > shared->nzc[ 0 ][ h2 ] ) {
							const ULI temp = h2;
							Matrix2HilbertCoordinates::IntegerToHilbert( beta[ h2-start ][ shared->nzc[ 0 ][ h2 ] ].i() / max_m,
								beta[ h2-start ][ shared->nzc[ 0 ][ h2 ] ].j() / max_n,
								h1, h2 );
							assert( h1 == 0 );
							assert( h2 == temp );
							h2 = temp;
							assert( row / max_m == beta[ h2-start ][ shared->nzc[ 0 ][ h2 ] ].i() / max_m );
							assert( col / max_n == beta[ h2-start ][ shared->nzc[ 0 ][ h2 ] ].j() / max_n );
						}
#endif
						beta[ h2-start ][ --(shared->nzc[ 0 ][ h2 ]) ] = (*(shared->original))[ i ];
						assert( row == beta[ h2-start ][ shared->nzc[ 0 ][ h2 ] ].i() );
						assert( col == beta[ h2-start ][ shared->nzc[ 0 ][ h2 ] ].j() );
					}
				}

#ifndef NDEBUG
				std::cout << "Thread " << shared->id << ": " << smin_m << "," << smax_m << " times " << smin_n << "," << smax_n << std::endl;
#endif

				//sanity check
				for( size_t i=0; i<end-start; i++ ) {
					assert( shared->nzc[ 0 ][ i + start ] == 0 );
					if( beta[ i ].size() == 0 ) continue;
					if( ms[ i ] - minms[ i ] >= max_m ) {
						std::cerr << "BetaHilbert thread construction: rowwise range (" << (ms[ i ]) << " to " << minms[ i ] << ") over maximum size! (h2=" << (i+start) << ")" << std::endl;
						exit( 1 );
					}
					if( ns[ i ] - minns[ i ] >= max_n ) {
						std::cerr << "BetaHilbert thread construction: columnwise range over maximum size!" << std::endl;
						exit( 1 );
					}
				}
	
				//sync
				BetaHilbert::synchronise( mutex, cond, shared->sync, shared->P );
	
				//clear local temporary array
				delete [] nzc;
	
				//remove temporary values
				delete [] ms;
				delete [] ns;
				delete [] minms;
				delete [] minns;
			} else {
				std::cerr << "Invalid value for proceed: " << proceed << std::endl;
				exit( 1 );
			}

			//load into FBICRS
#ifndef NDEBUG
			std::cout << "Processor " << shared->id << ": loading into local FBICRS data structure..." << std::endl;
#endif
			FBICRS< T > dss( beta, m, n );
			beta.clear();

			//remember data use
			shared->bytes = dss.bytesUsed();

			//create local shadow of y to avoid write-contention
			T* y = NULL;
			if( BH_USE_GLOBAL_Y && id == 0 )
				std::cerr << "Warning: thread 0 will use global y as local y (see BH_USE_GLOBAL_Y)" << std::endl;
			if( !BH_USE_GLOBAL_Y || shared->id > 0 ) { //id 0 will take global output vector for local y
				y = new T[ shared->output_vector_size ];
				for( size_t i=0; i<shared->output_vector_size; i++ )
					y[ i ] = 0.0;
			}
	
			//make local y known to outside world
			shared->local_y = y;
	
			//exit construction mode
			shared->mode = 0;

			//signal end of construction
			pthread_mutex_lock( mutex );
			BetaHilbert::end( shared->end_mutex, shared->end_cond, shared->end_sync, shared->P );

			//enter daemon mode
#ifndef NDEBUG
			std::cout << "Thread " << id << " construction done, entering daemon mode." << std::endl;
#endif
			while( true ) {
				struct timespec clk_start, clk_stop; 
				pthread_cond_wait(  cond, mutex );
				pthread_mutex_unlock( mutex );

				if( shared->mode == 4 ) break;
	
				switch( shared->mode ) {
				case 5:
					if( !(BH_USE_GLOBAL_Y && shared->id == 0) ) {
						for( size_t i = 0; i < shared->output_vector_size; ++i ) {
							shared->local_y[ i ] = 0;
						}
					}
					break;
				case 3:
					assert( *(shared->input)  != NULL );
					assert( *(shared->output) != NULL );
					if( BH_USE_GLOBAL_Y && shared->id == 0 ) {
						y = *(shared->output);
						shared->local_y= y;
					}
					assert( y != NULL );

					clock_gettime( global_clock_id, &clk_start);
					shared->time = 0.0;
					for( size_t i=0; i<shared->repeat; i++ ) {
						dss.zxa_fb( *(shared->input), y );
						if( BH_COLLECT == -1 )
							break;
						else if( BH_COLLECT != 0 )
							BetaHilbert::synchronise( mutex, cond, shared->sync, shared->P );
						collectY( shared );
					}
					clock_gettime( global_clock_id, &clk_stop);
					shared->time  = (clk_stop.tv_sec-clk_start.tv_sec)*1000;
					shared->time += (clk_stop.tv_nsec-clk_start.tv_nsec)/1000000.0;
					break;
				case 2:
					assert( *(shared->input)  != NULL );
					assert( *(shared->output) != NULL );
					if( BH_USE_GLOBAL_Y && shared->id == 0 ) {
						y = (*shared->output);
						shared->local_y = y;
					}
					assert( y != NULL );
					clock_gettime( global_clock_id, &clk_start);
					shared->time = 0.0;
					for( size_t i=0; i<shared->repeat; i++ ) {
						dss.zax_fb( *(shared->input), y );
						if( BH_COLLECT == -1 )
							continue;
						else if( BH_COLLECT != 0 )
							BetaHilbert::synchronise( mutex, cond, shared->sync, shared->P );
						collectY( shared );
					}
					clock_gettime( global_clock_id, &clk_stop);
					shared->time  = (clk_stop.tv_sec-clk_start.tv_sec)*1000;
					shared->time += (clk_stop.tv_nsec-clk_start.tv_nsec)/1000000.0;
					break;
				default:
					std::cout << "Thread " << id << ": Error, undefined operation (" << shared->mode << ")!" << std::endl;
					exit( -1 );
				}
				shared->mode = 0;

				//signal end of operation
				pthread_mutex_lock( mutex );
				BetaHilbert::end( shared->end_mutex, shared->end_cond, shared->sync, shared->P );
			}

			//done
			if( !BH_USE_GLOBAL_Y || shared->id > 0 )
				delete [] y;
			return (NULL);
		}

		/**
		 * Code that collects a distributed output vector.
		 * 
		 * @param shared The local thread data.
		 */
		static void collectY( shared_data<T> *shared ) {
			const ULI s = shared->id;
			const ULI p = shared->P;
			const ULI m = shared->output_vector_size;
			shared_data<T> *datas = shared; datas -= s;

			if( BH_COLLECT == 0 ) {
				//do nothing, will be handled after SPMD code
			} else if( BH_COLLECT == 1 ) {
				//use block distribution
				const ULI blocksize = (m % p == 0 ) ? m / p : m / p + 1;
				const ULI m_start = s * blocksize;
				const ULI m_end   = m_start + blocksize > m ? m : m_start + blocksize;

				if( s == p-1 ) assert( m_end == shared->output_vector_size );
				//do collect items of own block
				for( size_t i = m_start; i<m_end; i++ ) {
					for( size_t k = BH_USE_GLOBAL_Y; k < p; k++ ) { //start at k=1 when local y of proc. 0 is global y itself
						(*(shared->output))[ i ] += datas[ k ].local_y[ i ];
					}
				}
			} else if( BH_COLLECT == 2 ) {
				//parallel reducing collect
				size_t step_p = BH_REDUCE_BASE;
				size_t prev_step_p = 1;
				const ULI blocksize = (m % p == 0 ) ? m / p : m / p + 1;
				const ULI m_start = s * blocksize;
				const ULI m_end   = m_start + blocksize > m ? m : m_start + blocksize;
				while( prev_step_p < p ) {
					for( size_t start_p = 0; start_p < p; start_p += step_p ) {
						//do collect items of own block
						for( size_t k = start_p + prev_step_p; k < p && k < start_p + step_p; k += prev_step_p ) {
							for( size_t i = m_start; i<m_end; i++ )
								datas[ start_p ].local_y[ i ] += datas[ k ].local_y[ i ];
						}
					}
					prev_step_p = step_p;
					step_p *= BH_REDUCE_BASE;
				}
				//if datas[0].local_y != BetaHilbert<T>::output, one final step is required:
				if( !BH_USE_GLOBAL_Y )
					for( size_t i = m_start; i<m_end; i++ )
						(*(shared->output))[ i ] += datas[ 0 ].local_y[ i ];
			} else if( BH_COLLECT == 3 ) {
				size_t step_p = 2;
				size_t prev_step_p = 1;
				while( prev_step_p < p ) {
					for( size_t start_p = 0; start_p < p; start_p += step_p ) {
						if( prev_step_p > 1 )
							BetaHilbert::synchronise( shared->mutex, shared->cond, shared->sync, shared->P );
						if( start_p != s ) continue;
						for( size_t k = start_p + prev_step_p; k < p && k < start_p + step_p; k += prev_step_p )
							for( size_t i = 0; i < m; i++ )
								datas[ start_p ].local_y[ i ] += datas[ k ].local_y[ i ];
					}
					prev_step_p = step_p;
					step_p *= 2;
				}
				//if datas[0].local_y != BetaHilbert<T>::output, one final step is required:
				if( !BH_USE_GLOBAL_Y && s == 0 )
					for( size_t i = 0; i < m; i++ )
						(*(shared->output))[ i ] += datas[ 0 ].local_y[ i ];
			} else {
				std::cerr << "Error: output vector collection strategy " << BH_COLLECT << " not implemented!" << std::endl;
				exit( 1 );
			}
		}

		/**
		 * Computes z=xA in place.
		 * 
		 * @param x The input vector of length m.
		 * @param z The output vector of length n.
		 */
		virtual void zxa( const T* x, T* z ) {
			zxa( x, z, 1 );
		}

		/**
		 * Computes z=xA in place, a given number of times successively.
		 * 
		 * @param x The input vector of length m.
		 * @param z The output vector of length n.
		 * @param repeat The number of times to repeat the in-place SpMV.
		 */
		virtual void zxa( const T* x, T* z, const size_t repeat ) {
			//set all daemon threads to do zxa
			for( size_t i=0; i<P; i++ ) {
				thread_data[ i ].mode = 3;
				thread_data[ i ].repeat = repeat;
			}

			//set input vector
			input = x;

			//set output vector
			output = z;

			//wake up all daemon threads
			pthread_mutex_lock( &end_mutex );
			pthread_mutex_lock( &mutex );
			pthread_cond_broadcast( &cond );
			pthread_mutex_unlock( &mutex );

			//wait for end of operation
			wait();

			//sequential final collect method
			if( BH_COLLECT == 0 )
				for( size_t i=0; i<this->nor; i++ )
					for( size_t s=BH_USE_GLOBAL_Y; s<P; s++ ) //start with s=1 when local y of proc 0 is global y (or rather z) itself
						z[ i ] += thread_data[ s ].local_y[ i ];

			//unset vectors
			input  = NULL;
			output = NULL;
		}

		/**
		 * Reset all local output vectors to zero.
		 */
		void reset() {
			//reset all local output vectors
			for( size_t i=0; i<P; i++ ) {
				thread_data[ i ].mode = 5;
			}

			//wake up all daemon threads
			pthread_mutex_lock( &end_mutex );
			pthread_mutex_lock( &mutex );
			pthread_cond_broadcast( &cond );
			pthread_mutex_unlock( &mutex );

			//wait for end of operation
			wait();
		}

		/**
		 * Overloaded mv call; allocates output vector using numa_interleaved.
		 * 
		 * @param x The input vector of length n.
		 * @return Newly allocated output vector of length m.
		 */
		virtual T* mv( const T* x ) {
			//over-allocate in case we use matrices that are too small to be
			//vectorised internally (assuming we are using vecBICRS as sub_ds).
			size_t allocsize = (this->nor + 1) * sizeof( T );
			//allocate, either using libnuma (to interleave for i86)
			//or using dynamic aligned allocs (for Xeon Phi MICs)
#ifndef _NO_LIBNUMA
			T* ret = (T*) numa_alloc_interleaved( allocsize );
#else
			T* ret = (T*) _mm_malloc( allocsize, 64 ); //instead of `new T[ nor ];'
#endif
			//set to 0-vector
			for( ULI i=0; i<this->nor; i++ ) {
				ret[ i ] = this->zero_element;
			}

			//reset output vectors
			reset();

			//do in-place SpMV
			zax( x, ret );

			//return new output vector
			return ret;
		}

		/**
		 * Computes z=Ax in place.
		 * 
		 * @param x The input vector of length n.
		 * @param z The output vector of length m.
		 */
		virtual void zax( const T* x, T* z ) {
			zax( x, z, 1, 0, NULL );
		}

		/**
		 * Computes z=Ax in place, a given number of times in succession, and measures the time taken.
		 * 
		 * @param x The input vector of length n.
		 * @param z The output vector of length m.
		 * @param repeat   The number of times this SpMV will be repeated.
		 * @param clock_id Which system clock to use for timing.
		 * @param elapsed_time Where to add the elapsed time to.
		 */
		virtual void zax( const T* x, T* z, const size_t repeat, const clockid_t clock_id, double *elapsed_time ) {
			//set all daemon threads to do zax
			for( size_t i=0; i<P; i++ ) {
				thread_data[ i ].mode = 2;
				thread_data[ i ].repeat = repeat;
			}

			//set global clock ID
			global_clock_id = clock_id;

			//set input vector
			input = x;

			//set output vector
			output = z;

			//wake up all daemon threads
			pthread_mutex_lock( &end_mutex );
			pthread_mutex_lock( &mutex );
			pthread_cond_broadcast( &cond );
			pthread_mutex_unlock( &mutex );

			//wait for end of operation
			wait();

			//sequential final collect method
			if( BH_COLLECT == 0 )
				for( size_t i=0; i<this->nor; i++ )
					for( size_t s=BH_USE_GLOBAL_Y; s<P; s++ ) //start with s=1 when local y of proc 0 is global y (or rather z) itself
						z[ i ] += thread_data[ s ].local_y[ i ];

			//get elapsed time
			double maxtime = 0.0;
			for( size_t i=0; i<P; i++ ) {
				const double curtime = thread_data[ i ].time;
				if( curtime > maxtime ) maxtime = curtime;
			}
			if( elapsed_time != NULL )
				*elapsed_time += maxtime;

			//unset vectors
			input  = NULL;
			output = NULL;
		}

		/** @see SparseMatrix::getFirstIndexPair */
		virtual void getFirstIndexPair( ULI &i, ULI &j ) {
			std::cerr << "Warning: BetaHilbert::getFirstIndexPair has no unique answer since it implements a parallel multiplication!\nIgnoring call..." << std::endl;
		}

		/** @return The number of bytes used by all threads corresponding to this scheme. */
		virtual size_t bytesUsed() {
			size_t ret = 0;
			for( size_t s = 0; s < P; ++s )
				ret += thread_data[ s ].bytes;
			return ret;
		}

};

template< typename T > size_t BetaHilbert< T >::P = 0;

template< typename T > clockid_t BetaHilbert< T >::global_clock_id = 0;

#endif

