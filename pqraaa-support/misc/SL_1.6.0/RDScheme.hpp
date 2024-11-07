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


/*
 * File created by:
 *     A. N. Yzelman, Dept. of Computer Science, KU Leuven, 2011.
 */


#include <iostream>
#include <vector>
#include <map>
#include <pthread.h>

#ifndef _NO_LIBNUMA
 #include <numa.h>
#endif

#include "SparseMatrix.hpp"
#include "Hilbert.hpp"

#ifndef _H_RDScheme
#define _H_RDScheme

/*
 * When defined, thread 0 will use the global y vector for its local
 * computations. This introduces extra work in the form of one sync,
 * and the amount of available processors for the parallel collect.
 * The advantage is less data replication of the output vector.
 *
 * The synchronisation cannot be prevented. Using all processors in
 * the collect code is possible but has not been programmed currently.
 */
//#define RDScheme_GLOBAL_Y

#ifndef _TESTMODE
 /*
  * When defined, RDScheme will not collect output results in the
  * global output vector passed to this library. Used for timing
  * the true SpMV speeds only.
  */
 #define RDScheme_NO_COLLECT
#endif

/** Shared data for RDScheme threads. */
template< typename T >
class RDScheme_shared_data {

	public:

		/** Thread ID */
		size_t id;

		/** Total number of processors */
		size_t P;

		/** 0 undef, 1 init, 2 zax, 3 zxa, 4 exit */
		unsigned char mode;

		/** how many times to repeat the operation set in `mode' (above, only for 2 and 3) */
		unsigned long int repeat;

		/** Array of vectors of thread-local nonzeroes. */
		std::vector< Triplet< T > > *original;

		/** Will store rowsums */
		size_t *nzb;

		/** Will store local timing */
		double time;

		/** Will store memory use. */
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

		/** Base constructor. Will initialise all to invalid values or NULL. */
		RDScheme_shared_data(): id( -1 ), P( -1 ), mode( 0 ), repeat( 0 ), original( NULL ), nzb( NULL ), time( 0 ), 
				mutex( NULL ), cond( NULL ), end_mutex( NULL ), end_cond( NULL ),
				sync( NULL ), end_sync( NULL ),
				output_vector_size( -1 ), output_vector_offset( -1 ) {}

		/**
		 * Default constructor.
		 *
		 * @param _id Thread ID.
		 * @param _P  Number of SPMD processes.
		 * @param _original Original set of nonzeroes, split by blocks.
		 * @param _nzb Number of sparse blocks.
		 * @param _mutex Sync mutex.
		 * @param _cond  Sync condition.
		 * @param _end_mutex End sync mutex.
		 * @param _end_cond   End sync condition.
		 * @param _sync Sync counter.
		 * @param _end_sync End sync counter.
		 * @param _ovsize Output vector size.
		 * @param _ovoffset Output vector start position (global view).
		 */
		RDScheme_shared_data( size_t _id, size_t _P,
				std::vector< Triplet< double > > *_original,
				size_t *_nzb, 
				pthread_mutex_t *_mutex, pthread_cond_t *_cond, pthread_mutex_t *_end_mutex, pthread_cond_t *_end_cond,
				size_t *_sync, size_t *_end_sync,
				size_t _ovsize, size_t _ovoffset ):
				id( _id ),  P( _P ), mode( 1 ), repeat( 1 ), original( _original ), nzb( _nzb ), time( 0 ),
				mutex( _mutex ), cond( _cond ), end_mutex( _end_mutex ), end_cond( _end_cond ),
				sync( _sync ), end_sync( _end_sync ),
				output_vector_size( _ovsize ), output_vector_offset( _ovoffset ), local_y( NULL ) {}
};

/** Full parallel row-distributed SpMV, based on CSB (Morton curve + Cilk) and PThreads.
 *  Inspired by Aydin & Gilbert's CSB, and comments by Patrick Amestoy on the BICRS Hilbert scheme. */
template< typename T, typename DS >
class RDScheme: public SparseMatrix< T, ULI > {

	private:

	protected:

		/** Number of threads to fire up */
		static size_t P;

		/** Input vector */
		static const T* input;

		/** Output vector */
		static T* output;

		/** p_threads associated to this data strcuture */
		pthread_t *threads;

		/** array of initial thread data */
		RDScheme_shared_data<T> *thread_data;

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

	public:

		/** Base constructor. Reads input from file. */
		RDScheme( const std::string file, T zero ) {
			this->loadFromFile( file, zero );
		}

		/** Base constructor. Reads input from a set of triplets. */
		RDScheme( std::vector< Triplet< T > >& input, ULI m, ULI n, T zero ) {
			load( input, m, n, zero );
		}

		/** Base deconstructor. */
		virtual ~RDScheme() {
			//set all daemon threads to exit mode
			for( size_t i=0; i<P; i++ )
				thread_data[ i ].mode = 4;

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

		/**
		 * Lets the calling thread wait for the end of the SpMV multiply.
		 */
		void wait() {
			//wait for end signal
			pthread_cond_wait( &end_cond, &end_mutex );
			pthread_mutex_unlock( &end_mutex );
		}

		/**
		 * Loads a sparse matrix from an input set of triplets.
		 *
		 * @param input The input set of triplets.
		 * @param m     The number of rows.
		 * @param n	The number of columns.
		 * @param zero  What constitutes a zero in this sparse matrix instance.
		 */
		virtual void load( std::vector< Triplet< T > >& input, const ULI m, const ULI n, const T zero ) {
			//get number of cores available
			P = MachineInfo::getInstance().cores();

#ifndef _NO_LIBNUMA
			//set kernel to local thread allocation if it wasn't already the case
			numa_set_localalloc();
#endif

			//base settings
			this->zero_element = zero;
			this->nor = m;
			this->noc = n;
			this->nnz = input.size();

			size_t *nzb = new size_t [ this->m() ];
			for( size_t i=0; i<m; i++ ) nzb[ i ] = 0;

			//create P threads :)
			this->threads = new pthread_t[ P ];
			//initialize local initialisation data
			thread_data = new RDScheme_shared_data<T>[ P ];
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
				thread_data[ i ] = RDScheme_shared_data<T>( i, P, &input, nzb, &mutex, &cond, &end_mutex, &end_cond, &sync, &end_sync, -1, -1 );
				//set fixed affinity for threads
				cpu_set_t mask;
				CPU_ZERO( &mask );
				CPU_SET ( i, &mask );

				//TODO: use hwloc for better numa-aware pinning
				/*hwloc_topology_t topology;
				hwloc_topology_init ( &topology );
				hwloc_topology_load( topology );
				hwloc_bitmap_t cpuset;*/

				//prepare attributes
				pthread_attr_t attr;
				pthread_attr_init( &attr );
				//set fixed affinity in attribute, so that it starts binded immediately
				pthread_attr_setaffinity_np( &attr, sizeof( cpu_set_t ), &mask );
				//fire up thread
				pthread_create( &threads[i], &attr, &RDScheme::thread, (void*) &thread_data[i] );
				//free attr
				pthread_attr_destroy( &attr );
			}

			//wait for threads to finish initialisation
			wait();

			//delete temporary array
			delete [] nzb;
		}

		/** End synchronisation code. */
		static void end( pthread_mutex_t* mutex, pthread_cond_t* cond, size_t *sync, const size_t P ) {
			pthread_mutex_lock( mutex );
			(*sync)++;
			if( *sync == P ) {
				//only one thread is waiting on this condition, use signal
				pthread_cond_signal( cond );
				*sync = 0;
			}
			pthread_mutex_unlock( mutex );
		}

		/** Synchronises all threads. */
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

		/** SPMD code for each thread involved with parallel SpMV multiplication. */
		static void* thread( void *data ) {
			//get short-hand notation
			RDScheme_shared_data<T>* shared  = (RDScheme_shared_data<T>*)data;
			const size_t id  = shared->id;
			const size_t P   = shared->P;
			const size_t nnz = shared->original->size();
			pthread_mutex_t *mutex      = shared->mutex;
			pthread_cond_t  *cond       = shared->cond;

			cpu_set_t mask;
			CPU_ZERO( &mask );
			pthread_getaffinity_np( pthread_self(), sizeof( cpu_set_t ), &mask );

			//sanity checks
			if( !CPU_ISSET( id, &mask ) ) {
				std::cerr << "Incorrect pinning for thread " << id << "!" << std::endl;
				exit( 1 );
			}
			for( size_t s=0; s<P; s++ ) {
				if( s==id ) continue;
				if( CPU_ISSET( s, &mask ) ) {
					std::cerr << "Thread " << id << " mask is larger than one core" << " (" << s << " is set)!" << std::endl;
					exit( 1 );
				}
			}

			//prepare to get global matrix dimensions
			ULI m, n;
			m = n = 0;
			//put rowsums in nzb
			const size_t blocksize = (nnz % P) > 0 ? nnz / P + 1 : nnz / P;
			for( size_t i=0; i<nnz; i++ ) {
				const unsigned long int currow = (*(shared->original))[ i ].i();
				const unsigned long int curcol = (*(shared->original))[ i ].j();
				if( currow >= id * blocksize && currow < (id + 1) * blocksize )
					shared->nzb[ currow ]++;
				if( currow > m ) m = currow;
				if( curcol > n ) n = curcol;
			}
			
			//dimensions are one higher than max indices
			m++;
			n++;

			//sync
			RDScheme::synchronise( mutex, cond, shared->sync, shared->P );

			//determine distribution
			const size_t nnz_target = nnz / P;
			size_t cursum = 0;

			//first sanity check
			for( unsigned long int i=0; i<m; i++ ) cursum += shared->nzb[ i ];
			assert( cursum == nnz );
			
			//continue
			cursum = 0;
			size_t start, end, k = 0;
			start = end = -1;
			//get start position for s=0 correct
			if( id == 0 ) start = 0;
			//do greedy load balancing to get ranges for prcoessors 0 to P-1
			for( size_t i = 0; i < m; i++ ) {
				cursum += shared->nzb[ i ];	
				if( cursum >= nnz_target ) {
					if( k == id ) end   = i + 1;
					if(k+1== id ) start = i + 1;
					k++;
					cursum = 0;
				}
			}
			//see if we missed out on any nonzeroes, and set an empty range if so
			if( start == static_cast< size_t >(-1) ) start = m;
			if(  end  == static_cast< size_t >(-1) ) end   = m; 
			//get end position for s=P-1 correct
			if( id == P-1 ) end = m;
			//derive output vector sizes
			shared->output_vector_size   = end - start;
			shared->output_vector_offset = start;
			assert( shared->output_vector_size <= m );
			assert( shared->output_vector_offset + shared->output_vector_size <= m );

			//copy to local first			
			std::vector< Triplet< T > > local;
			for( size_t i = 0; i < static_cast< size_t >(nnz); i++ ) {
				const size_t currow = (*(shared->original))[ i ].i();
				if( currow >= start && currow < end )
					local.push_back(
						Triplet< T >( (*(shared->original))[ i ].i() - start,
							(*(shared->original))[ i ].j(),
							(*(shared->original))[ i ].value )
					);
			}
			m = shared->output_vector_size; //new matrix size is new m times old n

			//load into datastructure
			DS dss( local, m, n, 0 );

			//remember memory usage
			shared->bytes = dss.bytesUsed();

			//create local shadow of y to avoid write-contention
			T* y = NULL;
#ifdef RDScheme_GLOBAL_Y
			if( id > 0 ) {
#endif
				if( shared->output_vector_size > 0 ) {
					y = new T[ shared->output_vector_size ];
					for( size_t i=0; i<shared->output_vector_size; i++ )
						y[ i ] = 0.0;
				} else
					y = NULL;
#ifdef RDScheme_GLOBAL_Y
			}
#endif
			shared->local_y = y;
	
			//exit construction mode
			shared->mode = 0;

			//signal end of construction
			pthread_mutex_lock( mutex );
			RDScheme::end( shared->end_mutex, shared->end_cond, shared->end_sync, shared->P );

			//enter daemon mode
			while( true ) {
				struct timespec clk_start, clk_stop; 
				pthread_cond_wait(  cond, mutex );
				pthread_mutex_unlock( mutex );

				if( shared->mode == 4 ) break;
	
#ifndef NDEBUG
				const double * const p_input  = RDScheme<T,DS>::input;
				const double * const p_output = RDScheme<T,DS>::output;
#endif
				switch( shared->mode ) {
				case 3:
					assert( p_input  != NULL );
					assert( p_output != NULL );
#ifdef RDScheme_GLOBAL_Y
					if( id == 0 ) {
						y = RDScheme::output;
						shared->local_y = y;
					}
#endif
					assert( y != NULL );

					clock_gettime( global_clock_id, &clk_start);
					shared->time = 0.0;
					for( unsigned long int i=0; i<shared->repeat; ++i )
						dss.zxa( RDScheme<T,DS>::input, y );
					clock_gettime( global_clock_id, &clk_stop);
					shared->time  = (clk_stop.tv_sec-clk_start.tv_sec)*1000;
					shared->time += (clk_stop.tv_nsec-clk_start.tv_nsec)/1000000.0;

#ifndef RDScheme_NO_COLLECT
					collectY( shared );
#endif
					break;
				case 2:
					assert( p_input  != NULL );
					assert( p_output != NULL );
#ifdef RDScheme_GLOBAL_Y
					if( id == 0 ) {
						y = RDScheme::output;
						shared->local_y = y;
					}
#endif
					assert( y != NULL );

					clock_gettime( global_clock_id, &clk_start);
					shared->time = 0.0;
					for( unsigned long int i=0; i<shared->repeat; ++i )
						dss.zax( RDScheme<T,DS>::input, y );
					clock_gettime( global_clock_id, &clk_stop);
					shared->time  = (clk_stop.tv_sec-clk_start.tv_sec)*1000;
					shared->time += (clk_stop.tv_nsec-clk_start.tv_nsec)/1000000.0;

#ifndef RDScheme_NO_COLLECT
					collectY( shared );
#endif
					break;
				default:
					std::cout << "Thread " << id << ": Error, undefined operation (" << shared->mode << ")!" << std::endl;
					exit( -1 );
				}
				shared->mode = 0;

				//signal end of operation
				pthread_mutex_lock( mutex );
				RDScheme::end( shared->end_mutex, shared->end_cond, shared->sync, shared->P );
			}

			//done
#ifdef RDScheme_GLOBAL_Y
			if( id != 0 )
#endif
				delete [] y;
			return (NULL);
		}

		/**
		 * Reduces a distributed output vector set into a single
		 * contiguous output vector at process 0.
		 *
		 * @param shared Which set of output vectors to reduce.
		 */
		static void collectY( RDScheme_shared_data<T> *shared ) {

#ifdef RDScheme_GLOBAL_Y
			//FIXME It could be possible to distribute work over all processors
			//instead of p-1 processors, but this requires some extra balancing.
			const size_t s = shared->id;
			if( s == 0 ) return;
#endif

			//do collect items of own block
			for( size_t i = 0; i < shared->output_vector_size; i++ ) {
#ifndef NDEBUG
				const double * const p_output = RDScheme<T,DS>::output;
				assert( p_output != NULL );
				assert( shared->local_y != NULL );
#endif
				RDScheme<T,DS>::output[ shared->output_vector_offset + i ] += shared->local_y[ i ];
			}
		}

#ifndef _NO_LIBNUMA
		/** Overloaded mv call; allocates output vector using numa_interleaved. */
		virtual T* mv( const T* x ) {
			T* ret = (T*) numa_alloc_interleaved( this->nor * sizeof( T ) );
			for( ULI i=0; i<this->nor; i++ ) ret[ i ] = this->zero_element;
			zax( x, ret );
			return ret;
		}
#endif

		/** @see SparseMatrix::zxa */
		virtual void zxa( const T* x, T* z ) {
			zxa( x, z, 1 );
		}

		/** @see SparseMatrix::zxa */
		virtual void zxa( const T* x, T* z, const unsigned long int repeat ) {
			//set all daemon threads to do zxa
			for( size_t i=0; i<P; i++ ) {
				thread_data[ i ].mode   = 3;
				thread_data[ i ].repeat = repeat;
			}

			//set input vector
			RDScheme<T,DS>::input = x;

			//set output vector
			RDScheme<T,DS>::output = z;

			//wake up all daemon threads
			pthread_mutex_lock( &end_mutex );
			pthread_mutex_lock( &mutex );
			pthread_cond_broadcast( &cond );
			pthread_mutex_unlock( &mutex );

			//wait for end of operation
			wait();

			//unset vectors
			RDScheme<T,DS>::input  = NULL;
			RDScheme<T,DS>::output = NULL;
		}

		/** See SparseMatrix::zax */
		virtual void zax( const T* x, T* z ) {
			zax( x, z, 1, 0, NULL );
		}

		/** See SparseMatrix::zax */
		virtual void zax( const T* x, T* z, const unsigned long int repeat, const clockid_t clock_id, double *elapsed_time ) {
			//set all daemon threads to do zax
			for( size_t i=0; i<P; i++ ) {
				thread_data[ i ].mode   = 2;
				thread_data[ i ].repeat = repeat;
			}

			//set global clock ID
			global_clock_id = clock_id;

			//set input vector
			RDScheme<T,DS>::input = x;

			//set output vector
			RDScheme<T,DS>::output = z;

			//wake up all daemon threads
			pthread_mutex_lock( &end_mutex );
			pthread_mutex_lock( &mutex );
			pthread_cond_broadcast( &cond );
			pthread_mutex_unlock( &mutex );

			//wait for end of operation
			wait();

			//get elapsed time
			double maxtime = 0.0;
			for( size_t i=0; i<P; i++ ) {
				const double curtime = thread_data[ i ].time;
				if( curtime > maxtime ) maxtime = curtime;
			}
			if( elapsed_time != NULL )
				*elapsed_time += maxtime;

			//unset vectors
			RDScheme<T,DS>::input  = NULL;
			RDScheme<T,DS>::output = NULL;
		}

		/** @return The memory usage of this data structure (summed over all threads). */
		virtual size_t bytesUsed() {
			size_t ret = 0;
			for( size_t s = 0; s < P; ++s )
				ret += thread_data[ s ].bytes;
			return ret;
		}

		/**
		 * Function disabled for parallel schemes!
		 * @see SparseMatrix::getFirstIndexPair
		 */
		virtual void getFirstIndexPair( ULI &i, ULI &j ) {
			std::cerr << "Warning: RDScheme::getFirstIndexPair has no unique answer since it implements a parallel multiplication!\nIgnoring call..." << std::endl;
		}
};

template< typename T, typename DS > size_t RDScheme< T, DS >::P = 0;

template< typename T, typename DS > const T* RDScheme< T, DS >::input  = NULL;

template< typename T, typename DS > T* RDScheme< T, DS >::output = NULL;

template< typename T, typename DS > clockid_t RDScheme< T, DS >::global_clock_id = 0;

#endif

