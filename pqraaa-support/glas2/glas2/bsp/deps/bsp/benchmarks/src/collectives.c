
#ifndef PRIMITIVE
/**
 * Controls which primitive is begin benchmarked.
 * 
 * 0: bsp_put (default)
 * 1: bsp_hpput
 * 2: bsp_get
 * 3: bsp_hpget
 * 10: MPI_Put
 * 11: MPI_Get
 * 12: MPI collectives (MPI_Allgather, MPI_Gather, and MPI_Bcast)
 */
 #define PRIMITIVE 0
#endif

#if PRIMITIVE < 4
 #include <mcbsp.h>
#elif PRIMITIVE > 9
 #include <mpi.h>
#endif

#if PRIMITIVE > 9 && PRIMITIVE < 12
 #define MPI_PREFIX window,
#else
 #define MPI_PREFIX
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INITIAL_REP 10000

const size_t  n1 = 1ul<<14; //16kb
const size_t  n2 = 1ul<<17; //128kB
const size_t  n3 = 1ul<<22; //4MB
const size_t  n4 = 1ul<<27; //128MB

double time( void ) {
#if PRIMITIVE < 4
	return bsp_time();
#elif PRIMITIVE > 9
	return MPI_Wtime();
#endif
}

void sync( void ) {
#if PRIMITIVE < 4
	bsp_sync();
#elif PRIMITIVE > 9
	MPI_Barrier( MPI_COMM_WORLD );
#endif
}

double max_time( const double time, const size_t s, const size_t p ) {
#if PRIMITIVE < 4
	bsp_send( 0, NULL, &time, sizeof(double) );
	double maximum = time;
	bsp_sync();
	if( s == 0 ) {
		unsigned int messages;
		const double * curtime;
		const void * tagpointer;
		bsp_qsize( &messages, NULL );
		if( messages != p ) {
			bsp_abort( "Not enough BSMP messages!\n" );
		}
		for( size_t k = 0 ; k < p; ++k ) {
			bsp_hpmove( (const void**)(&tagpointer), (const void**)(&curtime) );
			if( *curtime > maximum ) {
				maximum = *curtime;
			}
		}
		return maximum;
	}
	return 0;
#elif PRIMITIVE > 9
	double maximum = 0;
	MPI_Reduce( &time, &maximum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
	return maximum;
#endif
}

//all-to-all
void kernel(
#if PRIMITIVE > 9 && PRIMITIVE < 12
	MPI_Win window,
#endif
	const void * const chunk1, void * const chunk2, const size_t rep, const size_t n, const size_t s, const size_t p ) {
	const double t0 = time();
	for( size_t i = 0; i < rep; ++i ) {
#if PRIMITIVE > 9 && PRIMITIVE < 12
		//ignoring return value; on error will abort MPI job
		MPI_Win_start( MPI_COMM_WORLD, 0, window );
#endif
#if PRIMITIVE == 12
		MPI_Allgather( chunk1, n, MPI_CHAR, chunk2, n, MPI_CHAR, MPI_COMM_WORLD );
#else
		for( size_t k = 0; k < p; ++k ) {
 #if PRIMITIVE == 0
			bsp_put( k, chunk1, chunk2, s * n, n );
 #elif PRIMITIVE == 1
			bsp_hpput( k, chunk1, chunk2, s * n, n );
 #elif PRIMITIVE == 2
			bsp_get( k, chunk1, 0, chunk2 + s*n, n );
 #elif PRIMITIVE == 3
			bsp_hpget( k, chunk1, 0, chunk2 + s*n, n );
 #elif PRIMITIVE == 10
			//ignoring return value; on error will abort MPI job
			MPI_Put( chunk1, n, MPI_CHAR, (int)k, s * n, n, MPI_CHAR, window );
 #elif PRIMITIVE == 11
			MPI_Get( chunk2 + s*n, n, MPI_CHAR, (int)k, 0, n, MPI_CHAR, window );
 #endif
		}
 #if PRIMITIVE < 4 
		bsp_sync();
 #elif PRIMITIVE > 9
		//ignoring return value; on error will abort MPI job
		MPI_Win_complete( window );
 #endif
#endif
	}
	const double local_elapsed_time = time() - t0;
	const double global_elapsed_time = 1000.0 * max_time( local_elapsed_time, s, p ) / ((double)rep); //avg time in ms
	if( s == 0 ) {
		printf( "Average time taken: %lf ms. Experiments were repeated %ld times.\n\n", global_elapsed_time, rep );
	}
}

//all-to-one
void kernel2(
#if PRIMITIVE > 9 && PRIMITIVE < 12
	MPI_Win window,
#endif
	const void * const chunk1, void * const chunk2, const size_t rep, const size_t n, const size_t s, const size_t p ) {
	const double t0 = time();
	for( size_t i = 0; i < rep; ++i ) {
#if PRIMITIVE == 0
		bsp_put( 0, chunk1, chunk2, s * n, n );
		bsp_sync();
#elif PRIMITIVE == 1
		bsp_hpput( 0, chunk1, chunk2, s * n, n );
		bsp_sync();
#elif PRIMITIVE == 2
		for( size_t k = 0; s == 0 && k < p; ++k ) {
			bsp_get( k, chunk1, 0, chunk2 + s*n, n );
		}
		bsp_sync();
#elif PRIMITIVE == 3
		for( size_t k = 0; s == 0 && k < p; ++k ) {
			bsp_hpget( k, chunk1, 0, chunk2 + s*n, n );
		}
		bsp_sync();
#elif PRIMITIVE == 10
		//ignoring return value; on error will abort MPI job
		MPI_Win_start( MPI_COMM_WORLD, 0, window );
		MPI_Put( chunk1, n, MPI_CHAR, 0, s * n, n, MPI_CHAR, window );
		MPI_Win_complete( window );
#elif PRIMITIVE == 11
		MPI_Win_start( MPI_COMM_WORLD, 0, window );
		for( size_t k = 0; s == 0 && k < p; ++k ) {
			MPI_Get( chunk2 + s*n, n, MPI_CHAR, (int)k, 0, n, MPI_CHAR, window );
		}
		MPI_Win_complete( window );
#elif PRIMITIVE == 12
		MPI_Gather( chunk1, n, MPI_CHAR, chunk2, n, MPI_CHAR, 0, MPI_COMM_WORLD );
#endif
	}
	const double local_elapsed_time = time() - t0;
	const double global_elapsed_time = 1000.0 * max_time( local_elapsed_time, s, p ) / ((double)rep); //avg time in ms
	if( s == 0 ) {
		printf( "Average time taken: %lf ms. Experiments were repeated %ld times.\n\n", global_elapsed_time, rep );
	}
}

//one-to-all
void kernel3(
#if PRIMITIVE > 9 && PRIMITIVE < 12
	MPI_Win window,
#endif
	const void * const chunk1, void * const chunk2, const size_t rep, const size_t n, const size_t s, const size_t p ) {
	const double t0 = time();
	for( size_t i = 0; i < rep; ++i ) {
#if PRIMITIVE > 9 && PRIMITIVE < 12
		//ignoring return value; on error will abort MPI job
		MPI_Win_start( MPI_COMM_WORLD, 0, window );
#endif

#if PRIMITIVE == 2
		bsp_get( 0, chunk1, 0, chunk2, n );
#elif PRIMITIVE == 3
		bsp_hpget( 0, chunk1, 0, chunk2, n );
#elif PRIMITIVE == 0
		for( size_t k = 0; s == 0 && k < p; ++k ) {
			bsp_put( k, chunk1, chunk2, 0, n );
		}
#elif PRIMITIVE == 1
		for( size_t k = 0; s == 0 && k < p; ++k ) {
			bsp_hpput( k, chunk1, chunk2, 0, n );
		}
#elif PRIMITIVE == 10
		for( size_t k = 0; s == 0 && k < p; ++k ) {
			//ignoring return value; on error will abort MPI job
			MPI_Put( chunk1, n, MPI_CHAR, (int)k, 0, n, MPI_CHAR, window );
		}
#elif PRIMITIVE == 11
		for( size_t k = 0; s == 0 && k < p; ++k ) {
			MPI_Get( chunk2, n, MPI_CHAR, 0, 0, n, MPI_CHAR, window );
		}
#elif PRIMITIVE == 12
		//expand some effort to make Bcast broadcast into a different memory region
		//including at the process root.
		if( s == 0 ) {
			memcpy( chunk2, chunk1, n );
		}
		MPI_Bcast( chunk2, n, MPI_CHAR, 0, MPI_COMM_WORLD );
#endif
#if PRIMITIVE < 4
		bsp_sync();
#elif PRIMITIVE > 9 && PRIMITIVE < 12
		//ignoring return value; on error will abort MPI job
		MPI_Win_complete( window );
#endif
	}
	const double local_elapsed_time = time() - t0;
	const double global_elapsed_time = 1000.0 * max_time( local_elapsed_time, s, p ) / ((double)rep); //avg time in ms
	if( s == 0 ) {
		printf( "Average time taken: %lf ms. Experiments were repeated %ld times.\n\n", global_elapsed_time, rep );
	}
}

int main( int argc, char** argv ) {
#if PRIMITIVE > 9
	MPI_Init( &argc, &argv );
	int p, s;
	MPI_Comm_rank( MPI_COMM_WORLD, &s );
	MPI_Comm_size( MPI_COMM_WORLD, &p );
 #if PRIMITIVE == 10
	const char method[] = "MPI_Put";
 #elif PRIMITIVE == 11
	const char method[] = "MPI_Get";
 #elif PRIMITIVE == 12
	const char method[] = "MPI collectives";
 #endif
#elif PRIMITIVE < 4
	bsp_begin( bsp_nprocs() );
	const size_t p = bsp_nprocs();
	const size_t s = bsp_pid();
 #if PRIMITIVE == 0
	const char method[] = "bsp_put";
 #elif PRIMITIVE == 1
	const char method[] = "bsp_hpput";
 #elif PRIMITIVE == 2
	const char method[] = "bsp_get";
 #elif PRIMITIVE == 3
	const char method[] = "bsp_hpget";
 #endif
#endif

	size_t rep = INITIAL_REP;

	if( s == 0 ) {
		printf( "This application will benchmark successive collective calls using %s. Results are averaged over multiple runs. Time taken for synchronised entry and exit of the communication phases are included in the timings. This SPMD section counts %ld BSP threads.\n\n", method, (unsigned long int)p );
	}

#if PRIMITIVE > 9
	void * chunk1, * chunk2;
	if( MPI_Alloc_mem( n4/p, MPI_INFO_NULL, &chunk1 ) || MPI_Alloc_mem( n4, MPI_INFO_NULL, &chunk2 ) ) {
		fprintf( stderr, "Could not allocate memory!\n" );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}
	MPI_Win window;
 #if PRIMITIVE == 10
	MPI_Win_create( chunk2, n4, 0, MPI_INFO_NULL, MPI_COMM_WORLD, &window );
 #elif PRIMITIVE == 11
	MPI_Win_create( chunk1, n4/p, 0, MPI_INFO_NULL, MPI_COMM_WORLD, &window );
 #endif
#elif PRIMITIVE < 4
	//get chunks of memory
	void * const chunk1 = malloc( n4 / p );
	void * const chunk2 = malloc(   n4   );
	if( chunk1 == NULL || chunk2 == NULL ) {
		bsp_abort( "Could not allocate memory!" );
	}
#endif
#if PRIMITIVE < 2
	bsp_push_reg( chunk2, n4 );
#elif PRIMITIVE < 4
	bsp_push_reg( chunk1, n4 / p );
#endif
	//L1 expirement
	if( s == 0 ) {
		printf( "All-to-all %s benchmark, %ld bytes of data per BSP thread (should hit L1).\n", method, n1 );
		fflush( stdout );
	}
	sync();
	kernel( MPI_PREFIX chunk1, chunk2, rep, n1/p, s, p );

	//L2 expirement
	if( s == 0 ) {
		printf( "All-to-all %s benchmark, %ld bytes of data per BSP thread (should hit L2).\n", method, n2 );
		fflush( stdout );
	}
	sync();
	kernel( MPI_PREFIX chunk1, chunk2, rep/10, n2/p, s, p );

	//L3 expirement
	if( s == 0 ) {
		printf( "All-to-all %s benchmark, %ld bytes of data per BSP thread (should hit L3).\n", method, n3 );
		fflush( stdout );
	}
	sync();
	kernel( MPI_PREFIX chunk1, chunk2, rep/100, n3/p, s, p );

	//RAM expirement
	if( s == 0 ) {
		printf( "All-to-all %s benchmark, %ld bytes of data per BSP thread (should hit RAM).\n", method, n4 );	
		fflush( stdout );
	}
	sync();
	kernel( MPI_PREFIX chunk1, chunk2, rep/1000, n4/p, s, p );

	//L1 expirement
	if( s == 0 ) {
		printf( "\nAll-to-one %s benchmark, %ld bytes of data per BSP thread (should hit L1).\n", method, n1/p );
		fflush( stdout );
	}
	rep *= 10;
	sync();
	kernel2( MPI_PREFIX chunk1, chunk2, rep, n1/p, s, p );

	//L2 expirement
	if( s == 0 ) {
		printf( "All-to-one %s benchmark, %ld bytes of data per BSP thread (should hit L2).\n", method, n2/p );
		fflush( stdout );
	}
	sync();
	kernel2( MPI_PREFIX chunk1, chunk2, rep/10, n2/p, s, p );

	//L3 expirement
	if( s == 0 ) {
		printf( "All-to-one %s benchmark, %ld bytes of data per BSP thread (should hit L3).\n", method, n3/p );
		fflush( stdout );
	}
	sync();
	kernel2( MPI_PREFIX chunk1, chunk2, rep/100, n3/p, s, p );

	//RAM expirement
	if( s == 0 ) {
		printf( "All-to-one %s benchmark, %ld bytes of data per BSP thread (should hit RAM).\n", method, n4/p );
		fflush( stdout );
	}
	sync();
	kernel2( MPI_PREFIX chunk1, chunk2, rep/1000, n4/p, s, p );

	//L1 expirement
	if( s == 0 ) {
		printf( "\nOne-to-all %s benchmark, %ld bytes of data per BSP thread (should hit L1).\n", method, n1/p );
		fflush( stdout );
	}
	sync();
	kernel3( MPI_PREFIX chunk1, chunk2, rep, n1/p, s, p );

	//L2 expirement
	if( s == 0 ) {
		printf( "One-to-all %s benchmark, %ld bytes of data per BSP thread (should hit L2).\n", method, n2/p );
		fflush( stdout );
	}
	sync();
	kernel3( MPI_PREFIX chunk1, chunk2, rep/10, n2/p, s, p );

	//L3 expirement
	if( s == 0 ) {
		printf( "One-to-all %s benchmark, %ld bytes of data per BSP thread (should hit L3).\n", method, n3/p );
		fflush( stdout );
	}
	sync();
	kernel3( MPI_PREFIX chunk1, chunk2, rep/100, n3/p, s, p );

	//RAM expirement
	if( s == 0 ) {
		printf( "One-to-all %s benchmark, %ld bytes of data per BSP thread (should hit RAM).\n", method, n4/p );
		fflush( stdout );
	}
	sync();
	kernel3( MPI_PREFIX chunk1, chunk2, rep/1000, n4/p, s, p );

	if( s== 0 ) {
		printf( "Cleaning up...\n" );
		fflush( stdout );
	}

#if PRIMITIVE < 4
	//free chunks
	free( chunk1 );
	free( chunk2 );
	bsp_end();
#elif PRIMITIVE > 9
 #if PRIMITIVE < 12
	MPI_Win_free( &window );
 #endif
	MPI_Free_mem( chunk1 );
	MPI_Free_mem( chunk2 );
	MPI_Finalize();
#endif
}

