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

#include<vector>
#include<iostream>
#include<cstdlib>

unsigned long int randuli( const unsigned long int m ) {
	return static_cast< unsigned long int >(  (static_cast< double >( rand() ) / static_cast< double >( RAND_MAX ) ) * static_cast< double >( m ) );
}

double randd() {
	return static_cast< double >( rand() ) / static_cast< double >( RAND_MAX );
}

int main() {
	unsigned long int m = 0;
	unsigned long int n = 0;
	double zero = 0.0;

	std::vector< Triplet< double > > naive = Triplet< double >::load( "test.trp", m, n );

	TS< double > test( naive, m, n, zero );
	
	double* x = new double[ n ];
	for( unsigned long int i=0; i<n; i++ )
		x[i]=1.0;
	std::cout << "Init done, starting MV..." << std::endl;
	double* z = test.mv( x );
	//std::cout << "[";
	double mean = 0;
	for( unsigned long int i=0; i<( m-1 ); i++ ) {
		mean += z[ i ] / m;
		//std::cout << z[ i ] << ",";
	}
	mean += z[ m-1 ] / m;
	//std::cout << z[ m-1 ] << "]" << std::endl;
	std::cout << "Mean: " << mean << std::endl;
	delete [] z;
	delete [] x;
}

