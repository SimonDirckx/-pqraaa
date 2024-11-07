
#include "mcbsp.h"
#include "stdio.h"

int main( int argc, char **argv ) {
	printf( "%u\n", bsp_nprocs() );
	return 0;
}

