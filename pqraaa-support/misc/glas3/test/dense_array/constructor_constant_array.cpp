#include <glas3/array/dense_array/type/constant_array.hpp>

#include <glas3/array/dense_array/algorithm/iostream.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j, k ;

	std::cout << "range construction" << '\n' ;
	std::cout << "glas3::constant_array<double> c1(-2.3, {2, 4, 3}) ;" << '\n' ;
	glas3::constant_array<double> c1(-2.3, {2, 4, 3}) ;
	std::cout << c1 << '\n' << '\n' ;
	if ( c1.size() != 24 || c1.ndof() != 0 || c1.shape().size() != 3 || c1.shape()[0] != 2 || c1.shape()[1] != 4 || c1.shape()[2] != 3 || c1.shape().shape().size() != 1 ) return 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 4; ++j ) {
    		for ( k = 0; k < 3; ++k ) {
    			if ( c1[i + 2 * j + 2 * 4 * k] != -2.3 || c1({i, j, k}) != -2.3 ) return 1 ;
    		}
    	}
    }

	std::cout << "range construction" << '\n' ;
	std::cout << "glas3::constant_array<double> c2(-2.3) ;" << '\n' ;
	glas3::constant_array<double> c2(-2.3) ;
	std::cout << c2 << '\n' << '\n' ;
	if ( c2.size() != 1 || c2.ndof() != 0 || c2.shape().size() != 0 || c2.shape().shape().size() != 1 ) return 1 ;
	if ( c2[0] != -2.3 ) return 1 ;

	std::cout << "range construction" << '\n' ;
	std::cout << "glas3::constant_array<double, 2> c3(-2.3, {2, 4}) ;" << '\n' ;
	glas3::constant_array<double, 2> c3(-2.3, {2, 4}) ;
	std::cout << c3 << '\n' << '\n' ;
	if ( c3.size() != 8 || c3.ndof() != 0 || c3.shape().size() != 2 || c3.shape()[0] != 2 || c3.shape()[1] != 4 || c3.shape().shape().size() != 1 ) return 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 4; ++j ) {
    		if ( c3[i + 2 * j] != -2.3 || c3({i, j}) != -2.3 ) return 1 ;
    	}
    }

	std::cout << "range construction" << '\n' ;
	std::cout << "glas3::constant_array<double, 1> c4(-2.3, {4}) ;" << '\n' ;
	glas3::constant_array<double, 1> c4(-2.3, {4}) ;
	std::cout << c4 << '\n' << '\n' ;
	if ( c4.size() != 4 || c4.ndof() != 0 || c4.shape().size() != 1 || c4.shape()[0] != 4 || c4.shape().shape().size() != 1 ) return 1 ;
	for ( i = 0; i < 4; ++i ) {
		if ( c4[i] != -2.3 || c4({i}) != -2.3 ) return 1 ;
	}

	std::cout << "range construction" << '\n' ;
	std::cout << "glas3::constant_array<double, 0> c5(-2.3) ;" << '\n' ;
	glas3::constant_array<double, 0> c5(-2.3) ;
	std::cout << c5 << '\n' << '\n' ;
	if ( c5.size() != 1 || c5.ndof() != 0 || c5.shape().size() != 0 || c5.shape().shape().size() != 1 ) return 1 ;
	if ( c5[0] != -2.3 ) return 1 ;

	std::cout << "range construction" << '\n' ;
	std::cout << "glas3::constant_array<double, 3> c6(-2.3, {2, 4, 3}) ;" << '\n' ;
	glas3::constant_array<double, 3> c6(-2.3, {2, 4, 3}) ;
	std::cout << c6 << '\n' << '\n' ;
	if ( c6.size() != 24 || c6.ndof() != 0 || c6.shape().size() != 3 || c6.shape()[0] != 2 || c6.shape()[1] != 4 || c6.shape()[2] != 3 || c6.shape().shape().size() != 1 ) return 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 4; ++j ) {
    		for ( k = 0; k < 3; ++k ) {
    			if ( c6[i + 2 * j + 2 * 4 * k] != -2.3 || c6({i, j, k}) != -2.3 ) return 1 ;
    		}
    	}
    }

	std::cout << '\n' ;
	return 0;
}
