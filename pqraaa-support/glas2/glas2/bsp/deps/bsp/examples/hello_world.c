
#include <mcbsp.h>
#include <stdio.h>

void spmd( void );

void spmd( void ) {
	bsp_begin( bsp_nprocs()-1 > 0 ? bsp_nprocs()-1 : 1 );
	printf( "Hello world from thread %d in SPMD section two!\n", bsp_pid() );
	bsp_end();
}

int main( int argc, char **argv ) {
	bsp_begin( bsp_nprocs() );
	printf( "Hello world from thread %d in SPMD section one!\n", bsp_pid() );
	bsp_end();
	
	printf( "Initialising thread is in-between parallel SPMD sections one and two...\n" );

	bsp_init( &spmd, argc, argv );
	spmd();

	printf( "Initialising thread is about to exit...\n" );
	return 0;
}

