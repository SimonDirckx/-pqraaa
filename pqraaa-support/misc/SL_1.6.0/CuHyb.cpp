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
 *     A. N. Yzelman, Dept. of Computer Science, KU Leuven, 2014.
 */


#include "CuHyb.hpp"

int CuHyb::compareTriplets( const void * left, const void * right ) {
	const Triplet< double > one = **( (Triplet< double > **)left );
	const Triplet< double > two = **( (Triplet< double > **)right );
	if ( one.j() < two.j() )
		return -1;
	if ( one.j() > two.j() )
		return 1;
	return 0;
}

CuHyb::CuHyb() {}

CuHyb::CuHyb( std::string file, double zero ) {
	this->loadFromFile( file, zero );
}
	
CuHyb::CuHyb( const ULI number_of_nonzeros, const ULI number_of_rows, const ULI number_of_cols, double zero ) {
	this->nnz = number_of_nonzeros;
	this->nor = number_of_rows;
	this->noc = number_of_cols;
	this->zero_element = zero;
}

CuHyb::CuHyb( std::vector< Triplet< double > > input, ULI m, ULI n, double zero ) {
	load( input, m, n, zero );
}

CuHyb::~CuHyb() {
	cudaFree( GPUx );
	cudaFree( GPUz );
	cusparseStatus_t status = cusparseDestroyHybMat( hybA );
	if( status != CUSPARSE_STATUS_SUCCESS ) {
		std::cerr << "Warning: error during destruction of CuSparse HYB matrix!" << std::endl;
	}
	status = cusparseDestroyMatDescr( descrA );
	if( status != CUSPARSE_STATUS_SUCCESS ) {
		std::cerr << "Warning: error during destruction of CuSparse matrix descriptor!" << std::endl;
	}
	status = cusparseDestroy( handle );
	if( status != CUSPARSE_STATUS_SUCCESS ) {
		std::cerr << "Warning: error during destruction of CuSparse handle!" << std::endl;
	}
}

void CuHyb::load( std::vector< Triplet< double > >& input, ULI m, ULI n, double zero ) {
	std::cout << "\tLoading in a vector of " << input.size() << " triplets into CuHyb..." << std::endl;

	this->zero_element = zero;
	//find nnz
	this->nnz = input.size();

	this->nor = m;
	this->noc = n;

	//sanity check
	assert( this->nnz > 0 );

	//build better datastructure
	std::vector< std::vector< Triplet< double >* > > ds_vec( this->nor );
	
	//move input there
	std::vector< Triplet< double > >::iterator in_it = input.begin();
	for( ; in_it != input.end(); ++in_it ) {
		//Triplet< double >* cur = &(*in_it);
		const ULI currow = in_it->i();
		const double value = in_it->value;
		if( value == this->zero_element ) { 
			this->nnz--;
			continue;
		}
		ds_vec.at( currow ).push_back( &(*in_it) );
	}

	//allocate arrays
	   int * row_start = new    int[ this->nor + 1 ];
	double * ds        = new double[ this->nnz ];
	   int * col_ind   = new    int[ this->nnz ];

	//make CRS
	ULI index = 0;
	for( ULI currow = 0; currow < this->nor; currow++ ) {
		row_start[ currow ] = index;
		if( ds_vec.at( currow ).size() == 0 ) continue;
		qsort( &( ds_vec.at( currow )[ 0 ] ), ds_vec.at( currow ).size(), sizeof( Triplet< double >* ), &compareTriplets );
		std::vector< Triplet< double >* >::iterator row_it = ds_vec.at( currow ).begin();
		for( ; row_it != ds_vec.at( currow ).end(); row_it++ ) {
			const Triplet< double > cur = *(*row_it);
			ds[ index ] = cur.value;
			col_ind[ index ] = cur.j();
			index++;
		}
	}
	row_start[ this->nor ] = this->nnz;
	//sanity check
	assert( index == this->nnz );

	//record top-left position
	i = 0;
	while( row_start[ i ] == row_start[ i+1 ] && i < this->nor ) {
		++i;
	}
	j = col_ind[ 0 ];

	//convert to HYB; use CuSparse library for this. First offload CRS structure
	double * GPUds;
	   int * GPUrow_start, * GPUcol_ind;
	cudaError_t cudaStat = cudaMalloc( (void**)&GPUrow_start, nnz * sizeof(int) );
	if( cudaStat != cudaSuccess ) {
		std::cerr << "Could not allocate GPU array!" << std::endl;
		exit( EXIT_FAILURE );
	}
	cudaStat = cudaMalloc( (void**)&GPUcol_ind, nnz * sizeof(int) );
	if( cudaStat != cudaSuccess ) {
		std::cerr << "Could not allocate GPU array!" << std::endl;
		exit( EXIT_FAILURE );
	}
	cudaStat = cudaMalloc( (void**)&GPUds, nnz * sizeof(double) ); 
	if( cudaStat != cudaSuccess ) {
		std::cerr << "Could not allocate GPU array!" << std::endl;
		exit( EXIT_FAILURE );
	}
	cudaStat = cudaMalloc( (void**)&GPUx, this->noc * sizeof(double) );
	if( cudaStat != cudaSuccess ) {
		std::cerr << "Could not allocate GPU array!" << std::endl;
		exit( EXIT_FAILURE );
	}
	cudaStat = cudaMalloc( (void**)&GPUz, this->nor * sizeof(double) );
	if( cudaStat != cudaSuccess ) {
		std::cerr << "Could not allocate GPU array!" << std::endl;
		exit( EXIT_FAILURE );
	}
	
	cudaMemcpy( GPUrow_start, row_start, (this->nor+1) * sizeof(int),    cudaMemcpyHostToDevice );
	cudaMemcpy( GPUcol_ind,   col_ind,   (this->nnz)   * sizeof(int),    cudaMemcpyHostToDevice );
	cudaMemcpy( GPUds,        ds,        (this->nnz)   * sizeof(double), cudaMemcpyHostToDevice );

	//initialise cusparse library
	cusparseStatus_t status = cusparseCreate( &handle );
	if( status != CUSPARSE_STATUS_SUCCESS ) {
		std::cerr << "CuSparse init failed!" << std::endl;
		exit( EXIT_FAILURE );
	}
	// create and setup matrix descriptor
	status = cusparseCreateMatDescr( &descrA ); 
	if ( status != CUSPARSE_STATUS_SUCCESS ) {
		std::cerr << "Could not create sparse matrix descriptor!" << std::endl;
		exit( EXIT_FAILURE );
	}
	cusparseSetMatType(      descrA, CUSPARSE_MATRIX_TYPE_GENERAL );
	cusparseSetMatDiagType(  descrA, CUSPARSE_DIAG_TYPE_NON_UNIT );
	cusparseSetMatIndexBase( descrA, CUSPARSE_INDEX_BASE_ZERO );  
	//note descrA.fillType is unset since unused
	cusparseHybPartition_t partitionType = CUSPARSE_HYB_PARTITION_AUTO;
	status = cusparseCreateHybMat(&hybA );
	if( status != CUSPARSE_STATUS_SUCCESS ) {
		std::cerr << "Could not create initial HybMat matrix type!" << std::endl;
		exit( EXIT_FAILURE );
	}
	cusparseDcsr2hyb( handle,
		static_cast< int >( this->nor ), static_cast< int >( this->noc ),
		descrA, GPUds, GPUrow_start, GPUcol_ind,
		hybA, 0, partitionType ); //0 parameter is unused
		//wait until finished
	cudaDeviceSynchronize();
	//clean up GPU
	cudaFree( GPUds );
	cudaFree( GPUrow_start );
	cudaFree( GPUcol_ind );
	
	//clean up host
	delete [] row_start;
	delete [] ds;
	delete [] col_ind;
	//done
	std::cout << "\t" << index << " nonzeroes loaded into CuHyb structure, and offloaded." << std::endl;
}

void CuHyb::zxa( const double *__restrict__ x, double *__restrict__ z ) {
	std::cerr << "CuHyb does not implement the zxa (yet), sorry!" << std::endl;
	exit( EXIT_FAILURE );
}

void CuHyb::zax( const double *__restrict__ x, double *__restrict__ z ) {
	//offload x, z
	cudaMemcpy( GPUz, z, (this->nor) * sizeof(double), cudaMemcpyHostToDevice );
	cudaMemcpy( GPUx, x, (this->noc) * sizeof(double), cudaMemcpyHostToDevice );
	//do compute
	const double zero = 0.0;
	const double  one = 1.0;
	const cusparseStatus_t status = cusparseDhybmv( handle,
					CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
					descrA, hybA, GPUx,
					&zero, GPUz );
	if ( status != CUSPARSE_STATUS_SUCCESS ) {
		std::cerr << "Error during CuSparse SpMV!" << std::endl;
		exit( EXIT_FAILURE );
	}
	//retrieve z
	cudaMemcpy( z, GPUz, (this->nor) * sizeof(double), cudaMemcpyDeviceToHost );
}

void CuHyb::zax( const double *__restrict__ x, double *__restrict__ z, const unsigned long int repeat, const clockid_t clock_id, double *elapsed_time ) {
	//offload x, z
	cudaMemcpy( GPUz, z, (this->nor) * sizeof(double), cudaMemcpyHostToDevice );
	cudaMemcpy( GPUx, x, (this->noc) * sizeof(double), cudaMemcpyHostToDevice );

	//init
	struct timespec start, stop; 
	const double zero = 0.0;
	const double  one = 1.0;

	//start timer
	if( elapsed_time != NULL ) {
		clock_gettime( clock_id, &start);
	}
	//start repeated kernel
	for( unsigned long int i = 0; i < repeat; i++ ) {
		const cusparseStatus_t status = cusparseDhybmv(
							handle,
							CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
							descrA, hybA, GPUx,
							&zero, GPUz
						);
		assert( status == CUSPARSE_STATUS_SUCCESS );
	}
	//stop timer
	if( elapsed_time != NULL ) {
		clock_gettime( clock_id, &stop);
		double time = (stop.tv_sec-start.tv_sec)*1000;
		time += (stop.tv_nsec-start.tv_nsec)/1000000.0;
		*elapsed_time += time;
	}

	//retrieve z
	cudaMemcpy( z, GPUz, (this->nor) * sizeof(double), cudaMemcpyDeviceToHost );
}

size_t CuHyb::bytesUsed() {
	//return GPU-side memory used
	size_t available, total;
	cudaMemGetInfo( &available, &total );
	const size_t used = total - available;
	assert( total >= available );
	return used;
}	

void CuHyb::getFirstIndexPair( size_t &i, size_t &j ) {
	i = this->i;
	j = this->j;
}

