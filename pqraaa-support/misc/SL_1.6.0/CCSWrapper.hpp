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
 *     A. N. Yzelman, Dept. of Mathematics, Utrecht University, 2010.
 */


#include <vector>

#include "Triplet.hpp"
#include "SparseMatrix.hpp"

#ifndef _H_CCSWRAPPER
#define _H_CCSWRAPPER

/**
 * Automatically transforms a row-major scheme into an column-major scheme.
 * Can wrap around any SparseMatrix type, and automatically switches input
 * row indices with input column indices, and switches the z=Ax operation
 * with z=xA, and vice versa.
 */
template< typename T, typename SparseMatrixType, typename IND >
class CCSWrapper : public SparseMatrix< T, IND > {


	protected:

		/** Pointer to the underlying data structure. */
		SparseMatrixType *ds;

	protected:

		/**
		 * Helper function that transposes an input matrix in Triplet format,
		 * in-place.
		 *
		 * @param input The input matrix.
		 */
		void transposeVector( std::vector< Triplet< T > > &input ) {
			typename std::vector< Triplet< T > >::iterator it = input.begin();
			for( ; it != input.end(); ++it )
				it->transpose();
		}

	public:

		/** Default constructor (initialises with invalid data). */
		CCSWrapper() : ds( NULL ) {}

		/** Base file-based constructor. */
		CCSWrapper( std::string file, T zero = 0 ) { loadFromFile( file, zero ); }

		/** Base empty matrix constructor (sets nnz, rows, columns only). */
		CCSWrapper( const IND nnz, const IND nor, const IND noc, T zero ) {
			ds = new SparseMatrixType( nnz, noc, nor, zero );
		}

		/** Base Triplet-based constructor. */
		CCSWrapper( std::vector< Triplet< T > > &input, IND m, IND n, T zero ) {
			load( input, m, n, zero );
		}

		/** Base destructor. */
		virtual ~CCSWrapper() {
			if ( ds != NULL ) delete ds;
		}

		/** Triplet-based loader; first transposes, then calls nested constructor. */
		virtual void load( std::vector< Triplet< T > > &input, IND m, IND n, T zero ) {
			transposeVector( input );
			ds = new SparseMatrixType( input, n, m, zero );
		}

		/** File-based loader; reads file, then passes to Triplet-based loader. */
		virtual void loadFromFile( const std::string file, const T zero = 0 ) {
			ULI m, n;
			std::vector< Triplet< T > > VT = FileToVT::parse( file, m, n );
			load( VT, m, n, zero );
		}

		/** Returns the number of matrix rows (taking into account transposition). */
		virtual ULI m() { return ds->n(); }

		/** Returns the number of matrix columns (taking into account transposition). */
		virtual ULI n() { return ds->m(); }

		/** Returns the number of nonzeroes. */
		virtual ULI nzs() { return ds->nzs(); }

		/** @see SparseMatrix::getFirstIndexPair */
		virtual void getFirstIndexPair( IND &row, IND &col ) { ds->getFirstIndexPair( col, row ); }
		
		/** @see Matrix::mv */
		virtual T* mv( const T* x ) {
			T * const ret = new T[ m() ];
			for( unsigned long int i = 0; i < m(); i ++ ) ret[ i ] = ds->zero_element;
			ds->zxa( x, ret );
			return ret;
		}

		/** @see Matrix::zax (takes into account transposition) */
		virtual void zax( const T*__restrict__ x, T*__restrict__ z ) {
			ds->zxa( x, z );
		}

		/** @see Matrix::zxa (takes into account transposition) */
		virtual void zxa( const T*__restrict__ x, T*__restrict__ z ) {
			ds->zax( x, z );
		}

		/** @see Matrix::bytesUsed */
		virtual size_t bytesUsed() {
			return ds->bytesUsed();
		}
};

#endif

