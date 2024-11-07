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


/** \file
 * File created by:
 *     A. N. Yzelman, Dept. of Mathematics, Utrecht University, 2009.
 */

/*! \mainpage Sparse Library
 *
 * This is the doxygen code documentation of the sparse library. The latest version can be found at
 *     http://albert-jan.yzelman.net/software#SL
 *
 * Contact info is available via the same link. This doxygen is not exhaustive.
 * A good starting point on usage info is SparseMatrix. Which implementing subclass
 * (i.e., which sparse matrix storage format) is most suitable for you depends on
 * your application and usage scenario. For sequential computations a good baseline
 * might be CBICRS, for parallel computations the RDBHilbert class performs best
 * overall.
 *
 * This code is copyright by A. N. Yzelman;
 * Dept. of Mathematics, Utrecht University, 2007-2011;
 * Dept. of Computer Science, KU Leuven, 2011-2014.
 *
 * This library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.<br>
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.<br>
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses/
 */

#ifndef _H_M
#define _H_M

#include <time.h>

/**
 *  Defines operations common to all matrices, which are implemented in this library.
 */
template< typename T >
class Matrix {
	private:
	
	protected:

	public:

		/** Base constructor. */
		Matrix() {}

		/** Base deconstructor. */
		virtual ~Matrix() {}

		/** @return The number of rows. */
		virtual unsigned long int m() = 0;

		/** @return The number of columns. */
		virtual unsigned long int n() = 0;

		/** @return The number of nonzeros. */
		virtual unsigned long int nzs() {
			return m() * n();
		}

		/**
		 * Calculates z=Ax (where A is this matrix).
		 *
		 * @param x The input vector.
		 * @return The output vector z.
		 * @see Matrix::zax.
		 */
		virtual T* mv( const T* x ) = 0;

		/**
		 *  In-place z=Ax function.
		 *
		 *  @param x The x vector to multiply current matrix with.
		 *  @param z The result vector. Must be pre-allocated and its elements should be initialised to zero.
		 */
		virtual void zax( const T*__restrict__ x, T*__restrict__ z ) = 0;

		/** 
		 * Wrapper function to call the zax kernel multiple times successively,
		 * while timing the duration of the operation. Essentially computes
		 * z=A^k. The timing functionalities are optional.
		 * 
		 * @param x The input vector.
		 * @param z The output vector.
		 * @param k How many times to repeat the z=Ax operation.
		 * @param clock_id The POSIX realtime clock ID (optional).
		 * @param elapsed_time To which double the time taken should be added
		 *                     (optional; set to NULL to disable timing).
		 * @see Matrix::zax
		 */
		virtual void zax( const T*__restrict__ x, T*__restrict__ z, const size_t k, const clockid_t clock_id = 0, double *elapsed_time = NULL ) {
			struct timespec start, stop; 
			if( elapsed_time != NULL ) {
				clock_gettime( clock_id, &start);
			}
			for( size_t i = 0; i < k; ++i ) {
				zax( x, z );
			}
			if( elapsed_time != NULL ) {
				clock_gettime( clock_id, &stop);
				double time = (stop.tv_sec-start.tv_sec)*1000;
				time += (stop.tv_nsec-start.tv_nsec)/1000000.0;
				*elapsed_time += time;
			}
		}

		/**
		 * In-place Z=AX function, where A is m x n, Z = m x k, and X is n x k.
		 * The default implementation is to perform k successive z=Ax operations.
		 * Efficient schemes should override this default implementation to achieve
		 * higher performance.
		 * 
		 * @tparam k How many linear maps to compute. This is a template parameter so that the
		 *           compiler can optimise better.
		 * @param X The k input vectors of length n each.
		 * @param Z The k output vectors of length m each; each is assumed to be pre-initialised.
		 * @see Matrix::zax
		 */
		template< size_t k >
		void ZaX( const T*__restrict__ const *__restrict__ const X, T *__restrict__ const *__restrict__ const Z ) {
			//loop k times
			for( size_t i = 0; i < k; ++i ) {
				//call regular z=Ax function
				zax( X[i], Z[i] );
			}
			//done
		}

		/**
		 *  In-place z=xA function.
		 *
		 *  @param x The input vector to left-multiply to the current matrix.
		 *  @param z The output vector. Must be pre-allocated and initialised to zero for correct results.
		 */
		virtual void zxa( const T*__restrict__ x, T*__restrict__ z ) = 0;

		/**
		 * In-place Z=XA function, where A is m x n, Z = k x n, and X is k x m.
		 *
		 * Transposed variant of Matrix::ZaX, analogue to Matrix::zxa
		 * being the transposed operation of Matrix::zax.
		 * 
		 * @tparam k How many linear maps to compute. This is a template parameter so that the
		 *           compiler can optimise better.
		 * @param X The k input vectors of length m each.
		 * @param Z The k output vectors of length n each.
		 * @see Matrix::ZaX
		 * @see Matrix::zxa
		 */
		template< size_t k >
		void ZXa( const T *__restrict__ const *__restrict__ const X, T *__restrict__ const *__restrict__ const Z ) {
			//loop k times
			for( size_t i = 0; i < k; ++i ) {
				//call regular z=xA function
				zxa( X[i], Z[i] );
			}
			//done
		}

		/** Wrapper function to call the zxa kernel multiple times successively, while timing the operation duration. */
		virtual void zxa( const T*__restrict__ x, T*__restrict__ z, const unsigned long int repeat, const clockid_t clock_id = 0, double *elapsed_time = NULL ) {
			struct timespec start, stop; 
			if( elapsed_time != NULL ) clock_gettime( clock_id, &start);
			for( unsigned long int i=0; i<repeat; i++ )
				zxa( x, z );
			if( elapsed_time != NULL ) {
				clock_gettime( clock_id, &stop);
				double time = (stop.tv_sec-start.tv_sec)*1000;
				time += (stop.tv_nsec-start.tv_nsec)/1000000.0;
				*elapsed_time += time;
			}
		}

		/**
		 *  Function to query the amount of storage required by this sparse matrix.
		 *
		 *  @return The size of the sparse matrix in bytes.
		 */
		virtual size_t bytesUsed() = 0;
};

#endif //_H_M

