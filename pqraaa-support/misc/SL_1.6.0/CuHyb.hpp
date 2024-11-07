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


#include <assert.h>

#include <cuda_runtime.h>
#include <cusparse_v2.h>

#include "SparseMatrix.hpp"

#ifndef _H_CUHYB
#define _H_CUHYB

/** Wrapper class for the CuSparse HYB data structure for CUDA C. */
class CuHyb : public SparseMatrix< double, size_t > {

  private:

  protected:

	/** Handle to CuSparse. */
	cusparseHandle_t handle;

	/** CuSparse matrix descriptor. */
	cusparseMatDescr_t descrA;

	/** GPU-local matrix. */
	cusparseHybMat_t hybA;

	/** GPU-local buffer to the input vector. */
	double * GPUx;

	/** GPU-local buffer to the output vector. */
	double * GPUz;

 	/** Sorts 1D columnwise */
        static int compareTriplets( const void * left, const void * right );

	/** Top-left row coordinate. */
	size_t i;

	/** Top-right column coordinate. */
	size_t j;


  public:

	/** Base constructor. */
	CuHyb();

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	CuHyb( std::string file, double zero = 0 );
	
	/**
	 *  Base constructor which only initialises the internal arrays. Note that to gain a valid structure,
	 *  these arrays have to be filled by some external mechanism, and this mechanism needs to call CUDA
	 *  and CuSparse initialisation routines.
	 *  @param number_of_nonzeros The number of non-zeros to be stored.
	 *  @param number_of_rows The number of rows to be stored.
	 *  @param number_of_cols The number of columns of the matrix.
	 *  @param zero The element considered to be zero.
	 */
	CuHyb( const ULI number_of_nonzeros, const ULI number_of_rows, const ULI number_of_cols, double zero );

	/**
	 *  Constructor which transforms a collection of input triplets to CRS format.
	 *  The input collection is considered to have at most one triplet with unique
	 *  pairs of indeces. Unspecified behaviour occurs when this assumption is not
	 *  met.
	 *  @param input The input collection.
	 *  @param m Total number of rows.
	 *  @param n Total number of columns.
	 *  @param zero The element considered to be zero.
	 */
	CuHyb( std::vector< Triplet< double > > input, ULI m, ULI n, double zero );

	/** Base deconstructor. */
	virtual ~CuHyb();

	/** @see SparseMatrix::load */
	virtual void load( std::vector< Triplet< double > >& input, ULI m, ULI n, double zero );

	/** 
	 *  In-place z=xA function.
	 */
        virtual void zxa( const double *__restrict__ x, double *__restrict__ z );

	/** 
	 * In-place z=Ax function.
	 *  
	 * @param x The x vector to multiply current matrix with.
	 * @param z The result vector. Must be pre-allocated and its elements should be initialised to zero.
	 */
        virtual void zax( const double *__restrict__ x, double *__restrict__ z );

	/**
	 * In-place z=Ax function that supports repition and timing.
	 * 
	 * @param x The x vector to multiply current matrix with.
	 * @param z The result vector. Must be pre-allocated and its elements should be initialised to zero.
	 * @param repeat How many times to repeat the SpMV, essentially computing z=A^{repeat}x.
	 * @param clock_id Which POSIX realtime clock to use for timing.
	 * @param elapsed_time Where to add the elapsed time to.
	 */
	virtual void zax( const double *__restrict__ x, double *__restrict__ z, const unsigned long int repeat, const clockid_t clock_id = 0, double *elapsed_time = NULL );

	/** @see Matrix::bytesUsed */
	virtual size_t bytesUsed();

	/** @see SparseMatrix::getFirstIndexPair */
	virtual void getFirstIndexPair( size_t &i, size_t &j );

};

#endif

