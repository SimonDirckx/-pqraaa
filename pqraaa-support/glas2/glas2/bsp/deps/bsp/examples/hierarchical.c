
#include "mcbsp.h"

#include <stdio.h>
#include <stdlib.h>

void spmd1( void );
void spmd2_entry( void );
void spmd2_p0( double * const );

void spmd2_entry() {
	bsp_begin( bsp_nprocs()*2 );
	double local = 1.0;
	printf( "Thread %d has contribution %lf.\n", bsp_pid(), local );
	bsp_push_reg( &local, sizeof( double ) );
	bsp_sync();
	bsp_end();
}

void spmd2_p0( double * const a ) {
	bsp_begin( bsp_nprocs() );
	double *array = malloc( bsp_nprocs() * sizeof( double ) );
	bsp_push_reg( array,    bsp_nprocs() * sizeof( double ) );
	array[ 0 ] = 1.0;
	printf( "Thread 0 has contribution %lf.\n", array[ 0 ] );
	bsp_sync();
	printf( "Thread 0 combines partial results..." );
	for( unsigned int s = 0; s < bsp_nprocs(); ++s ) {
		bsp_direct_get( s, array, 0, array + s, sizeof( double ) );
	}
	for( unsigned int s = 0; s < bsp_nprocs(); ++s ) {
		*a += array[ s ];
	}
	printf( "Partial contribution: %lf.\n", *a );
	free( array );
	bsp_end();
}

void spmd1() {
	bsp_begin( 2 );
	double a = 0.0;
	bsp_push_reg( &a, sizeof( double ) );
	bsp_sync();
	printf( "Thread %d starts new BSP run...\n", bsp_pid() );
	bsp_init( &spmd2_entry, 0, NULL );
	spmd2_p0( &a );
	bsp_sync();
	if( bsp_pid() == 0 ) {
		printf( "Thread %d combines all results...", bsp_pid() );
		double temp;
		bsp_direct_get( 1, &a, 0, &temp, sizeof( double ) );
		a += temp;
		printf( " sum is %lf.\n", a );
	}
	bsp_end();
}

int main() {

	printf( "Sequential part 1\n" );

	bsp_init( &spmd1, 0, NULL );
	spmd1();

	printf( "Sequential part 2\n" );

	return EXIT_SUCCESS;	
}

