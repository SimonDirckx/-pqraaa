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


#include "TS.hpp"
#include "CRS.hpp"
#include "ICRS.hpp"
#include "Triplet.hpp"
#include "FileToVT.hpp"
#include "RDBHilbert.hpp"

#include <map>
#include <set>
#include <vector>
#include <sstream>
#include <iostream>

#define CLOCK_ID CLOCK_MONOTONIC

double elapsedTime( const struct timespec start, const struct timespec stop ) {
	double ret = (stop.tv_sec-start.tv_sec)*1000;
	ret += (stop.tv_nsec-start.tv_nsec)/1000000.0;
	return ret;
}

template< size_t k >
void verify( const size_t n, const double * const * const X, const double * const * const Z ) {
	//statistics
	double MSE, max_v, meanX, meanC;
	size_t max_i;
	//do verification for each vector separately
	for( size_t i = 0; i < k; ++i ) {
		//initialise
		MSE = meanX = meanC = 0.0;
		max_v = fabs(X[i][0] - Z[i][0]);
		max_i = 0;
		//loop over vector elements
		for( size_t j = 0; j < n; ++j ) {
			//current error and squared error
			const double cur_diff = (X[i][j] - Z[i][j]);
			const double cur_se   = cur_diff * cur_diff;
			//update means
			meanX += X[i][j] / static_cast< double >(n);
			meanC += Z[i][j] / static_cast< double >(n);
			//update MSE
			MSE += cur_se / static_cast< double >(n);
			//update maximum error
			if( fabs(cur_diff) > max_v ) {
				max_v = fabs(cur_diff);
				max_i = j;
			}
		}
		//report
		std::cout << "Verification on vector " << i << ": mean(X)=" << meanX << ", mean(C)=" << meanC << ", MSE=" << MSE << ", max error=" << max_v;
		std::cout << " when comparing X[" << max_i << "]=" << X[i][max_i] << " and C[" << max_i << "]=" << Z[i][max_i] << std::endl;
	}
	//done
}	

template< size_t k >
void resetVectors( const size_t m, const size_t n, double * const * const X, double * const * const Z ) {
	//use same seed for reproducability
	srand( 1337 );
	//initialise output vectors
	for( unsigned long int i = 0; i < m; ++i ) {
		for( size_t j = 0; j < k; ++j ) {
			Z[ j ][ i ] = 0.0;
		}
	}
	//initialise random input vectors
	for( unsigned long int j = 0; j < n; ++j ) {
		for( size_t i = 0; i < k; ++i ) {
			X[ i ][ j ] = 1.0 - rand() / static_cast< double >(RAND_MAX);
		}
	}
	//done
}

template< size_t k, class Matrix >
void verify_bench( const size_t inner_rep, Matrix &A, double **X, const double * const * const C, double **Z ) {
	//local data fields
	struct timespec start, stop; 
	double time;

	//initialise
	resetVectors< k >( A.m(), A.n(), X, Z );

	//do verify
	A.template ZaX< k >( X, Z );
	verify< k >( A.m(), Z, C );

	//if square matrix, do verification of ZXa as well
	if( A.m() == A.n() ) {
		std::cout << "\nExtra verification of Z=XA (only performed if A is square):\n";
		//new output vector alloc
		double **C2 = new double*[k];
		//for each of the k SpMVs
		for( size_t i = 0; i < k; ++i ) {
			//alloc new output
			C2[ i ] = new double[ A.n() ];
		}
		//init
		resetVectors< k >( A.n(), A.m(), X, C2 );
		resetVectors< k >( A.n(), A.m(), X, Z  );
		//do ZXa
		for( size_t i = 0; i < k; ++i ) {
			//delegate to z=xA
			A.zxa( X[i], C2[i] );
		}
		//specialised routine
		A.template ZXa< k >( X, Z );
		//call verification
		verify< k >( A.n(), Z, C2 );
		//cleanup
		for( size_t i = 0; i < k; ++i ) {
			delete [] C2[ i ];
		}
		delete [] C2;
	}

	//real run
	std::cout << "\nTimed run of " << inner_rep << " Z=AX computations:\n";
	clock_gettime( CLOCK_ID, &start);
	for( size_t i = 0; i < inner_rep; ++i ) {
		A.template ZaX< k >( X, Z );
	}
	clock_gettime( CLOCK_ID, &stop);
	time  = elapsedTime( start, stop );
	time /= static_cast< double >(inner_rep);
	std::cout << "Average time taken for one " << k << "-RHS SpMV multiplication is " << time << " ms." << std::endl;
	clock_gettime( CLOCK_ID, &start);
	for( size_t i = 0; i < inner_rep; ++i ) {
		for( size_t s = 0; s < k; ++s ) {
			A.zax( X[s], Z[s] );
		}
	}
	clock_gettime( CLOCK_ID, &stop);
	time  = elapsedTime( start, stop );
	time /= static_cast< double >(inner_rep);
	std::cout << "Average time taken for " << k << " consecutive SpMV multiplications, on independent input and output, is " << time << " ms." << std::endl;

	//done
}

std::vector< std::vector< Triplet< double > > > buildBetaHilbert( std::vector< Triplet< double > > &triplets, const unsigned long int m, const unsigned long int n ) {
	//will store the final FBICRS input
	std::vector< std::vector< Triplet< double > > > beta;
	//nonzero to block ID map
	std::vector< size_t > block_nonzero_id;
	//block ID to beta ID map
	std::map< size_t, size_t > block_beta;
	//Hilbert ID to block ID
	std::map< size_t, size_t > hilbert_block_id;
	//which block IDs are used
	std::set< size_t > block_used;
	//other useful constants, including the maximum number of block rows and block columns
	size_t h1, h2, c = 0, num_row_blocks = m / FBICRS< double >::beta_m, num_col_blocks = n / FBICRS< double >::beta_n;
	if( m % FBICRS< double >::beta_m > 0 ) {
		++num_row_blocks;
	}
	if( n % FBICRS< double >::beta_n > 0 ) {
		++num_col_blocks;
	}
	//for each nonzero
	for( size_t i = 0; i < triplets.size(); ++i ) {
		//get the block ID
		const size_t row = triplets[ i ].i() / FBICRS< double >::beta_m;
		const size_t col = triplets[ i ].j() / FBICRS< double >::beta_n;
		//derive 1D block index
		const size_t blk = row * num_col_blocks + col;
		//record index for this nonzero
		block_nonzero_id.push_back( blk );
		//record this block is used
		block_used.insert( blk );
	}
	//for each block that is used, derive its Hilbert ID
	for( const size_t block_id : block_used ) {
		//derive block coordinate
		const size_t row = block_id / num_col_blocks;
		const size_t col = block_id % num_col_blocks;
		//derive Hilbert coordinate
		Matrix2HilbertCoordinates::IntegerToHilbert( row, col, h1, h2 );
		//should fit within a single machine word
		assert( h1 == 0 );
		//store ID
		hilbert_block_id[ h2 ] = block_id;
	}
	//in the order of the Hilbert IDs, create blockID to final vector position map
	for( const std::pair< size_t, size_t > hilbertID_blockID : hilbert_block_id ) {
		//record position
		block_beta[ hilbertID_blockID.second ] = c++;
	}
	//do build beta
	beta.resize( c );
	//for each nonzero
	for( size_t i = 0; i < triplets.size(); ++i ) {
		//get block ID
		beta[ block_beta[ block_nonzero_id[ i ] ] ].push_back( triplets[i] );
	}
	//done, return beta
	return beta;
}

int main( int argc, char** argv ) {
	//local fields
	struct timespec start, stop; 
	double time;
	const size_t k = 4;
	const size_t default_inner_rep = 10;

	//check input
	if( argc == 1 || argc > 3 ) {
		std::cout << "Usage: " << argv[0] << "<MatrixMarket file> <rep (optional)>\n\n";
		std::cout << "Performs " << k << "-RHS SpMVs using various sequential formats. Will also perform\n";
		std::cout << "simple timing of each multiple-RHS SpMV by averaging speed over rep instances." << std::endl;
		return EXIT_SUCCESS;
	}

	//get filename
	const std::string file = std::string( argv[ 1 ] );

	//get inner_rep
	size_t inner_rep = default_inner_rep;
	if( argc >= 3 ) {
		std::istringstream arg( argv[2] );
		if( !(arg >> inner_rep) ) {
			inner_rep = default_inner_rep;
		}
	}

	//parse A in TS format
	TS< double > A_TS( file, 0 );

	//parse A in Triplet format
	std::vector< Triplet< double > > triplets = FileToVT::parse( file );

	//allocation
	double **X = new double*[k];
	double **C = new double*[k];
	double **Z = new double*[k];
	for( size_t i = 0; i < k; ++i ) {
		X[i] = new double[ A_TS.n() ];
		C[i] = new double[ A_TS.m() ];
	 	Z[i] = new double[ A_TS.m() ];
	}

	//verification run
	std::cout << "Verification of Z=AX:\n";
	resetVectors< k >( A_TS.m(), A_TS.n(), X, C );
	clock_gettime( CLOCK_ID, &start);
	for( size_t i = 0; i < k; ++i ) {
		A_TS.zax( X[i], C[i] );
	}
	clock_gettime( CLOCK_ID, &stop);
	time = elapsedTime( start, stop );
	std::cout << "Verification run completed in " << time << " ms." << std::endl;

	//test TS format
	std::cout << "\n==========================================\n";
	std::cout <<   "Triplet Scheme verification and benchmark:\n";
	std::cout <<   "==========================================\n\n";
	verify_bench< k >( inner_rep, A_TS, X, C, Z );

	//test CRS format
	std::cout << "\n==================================================\n";
	std::cout <<   "Compressed Row Storage verification and benchmark:\n";
	std::cout <<   "==================================================\n\n";
	CRS< double > A_CRS( triplets, A_TS.m(), A_TS.n(), 0 );
	verify_bench< k >( inner_rep, A_CRS, X, C, Z );

	//test ICRS format
	std::cout << "\n==============================================================\n";
	std::cout <<   "Incremental Compressed Row Storage verification and benchmark:\n";
	std::cout <<   "==============================================================\n\n";
	ICRS< double > A_ICRS( triplets, A_TS.m(), A_TS.n(), 0 );
	verify_bench< k >( inner_rep, A_ICRS, X, C, Z );

	//test Beta Hilbert format
	std::cout << "\n============================================================\n";
	std::cout <<   "Sequential Beta Hilbert (FBICRS) verification and benchmark:\n";
	std::cout <<   "============================================================\n\n";
	std::vector< std::vector< Triplet< double > > > beta = buildBetaHilbert( triplets, A_TS.m(), A_TS.n() );
	FBICRS< double > A_FBICRS( beta, A_TS.m(), A_TS.n(), 0 );
	verify_bench< k >( inner_rep, A_FBICRS, X, C, Z );

	//clean up
	for( size_t i = 0; i < k; ++i ) {
		delete [] X[i];
		delete [] C[i];
		delete [] Z[i];
	}
	delete [] X;
	delete [] C;
	delete [] Z;

	//done
	std::cout << std::endl;
	return EXIT_SUCCESS;
}

