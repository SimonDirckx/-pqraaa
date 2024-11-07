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


#include "Matrix2HilbertCoordinates.hpp"
#include <iostream>
#include <sstream>

int main( int argc, char **argv ) {
	if( argc < 2 ) {
		std::cout << "Usage: " << argv[ 0 ] << " <n>\n";
		std::cout << "\n\tprints an n times n grid of Hilbert coordinates to stdout." << std::endl;
		return EXIT_SUCCESS;
	}

	char *check = NULL;
	size_t n = strtoul( argv[1], &check, 10 );
	if( *check != '\0' ) {
		std::cerr << "Error parsing command line argument; is input an unsigned integer?" << std::endl;
		return EXIT_FAILURE;
	}

	for( size_t i = 0; i < n; ++i ) {
		std::cout << "\n";
		for( size_t j = 0; j < n; ++j  ) {
			size_t h1, h2;
			Matrix2HilbertCoordinates::IntegerToHilbert( i, j, h1, h2 );
			if( h1 == 0 ) {
				std::cout << h2 << " \t";
			} else {
				std::cout << h1 << h2 << " \t";
			}
		}
	}
	std::cout << "\n" << std::endl;
	return EXIT_SUCCESS;
}

