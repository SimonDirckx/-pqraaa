#include <glas3/array/dense_array/type/eye.hpp>

#include <glas3/array/dense_array/algorithm/iostream.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j, k ;

	std::cout << "eye construction" << '\n' ;
	std::cout << "glas3::eye<double> e1(2, 3) ;" << '\n' ;
	glas3::eye<double> e1(2, 3) ;
	std::cout << e1 << '\n' << '\n' ;
	if ( e1.size() != 8 || e1.ndof() != 0 || e1.shape().size() != 3 || e1.shape()[0] != 2 || e1.shape()[1] != 2 || e1.shape()[2] != 2 || e1.shape().shape().size() != 1 ) return 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 2; ++j ) {
    		for ( k = 0; k < 2; ++k ) {
    			if ( i == j && j == k ) { if ( e1[i + 2 * j + 2 * 2 * k] != 1 || e1({i, j, k}) != 1 ) return 1 ; }
    			else { if ( e1[i + 2 * j + 2 * 2 * k] != 0 || e1({i, j, k}) != 0 ) return 1 ; }
    		}
    	}
    }

	std::cout << "eye construction" << '\n' ;
	std::cout << "glas3::eye<double> e2({2, 4, 3}) ;" << '\n' ;
	glas3::eye<double> e2({2, 4, 3}) ;
	std::cout << e2 << '\n' << '\n' ;
	if ( e2.size() != 24 || e2.ndof() != 0 || e2.shape().size() != 3 || e2.shape()[0] != 2 || e2.shape()[1] != 4 || e2.shape()[2] != 3 || e2.shape().shape().size() != 1 ) return 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 4; ++j ) {
    		for ( k = 0; k < 3; ++k ) {
    			if ( i == j && j == k ) { if ( e2[i + 2 * j + 2 * 4 * k] != 1 || e2({i, j, k}) != 1 ) return 1 ; }
    			else { if ( e2[i + 2 * j + 2 * 4 * k] != 0 || e2({i, j, k}) != 0 ) return 1 ; }
    		}
    	}
    }

	std::cout << "eye construction" << '\n' ;
	std::cout << "glas3::eye<double, 2> e3(2, 2) ;" << '\n' ;
	glas3::eye<double, 2> e3(2, 2) ;
	std::cout << e3 << '\n' << '\n' ;
	if ( e3.size() != 4 || e3.ndof() != 0 || e3.shape().size() != 2 || e3.shape()[0] != 2 || e3.shape()[1] != 2 || e3.shape().shape().size() != 1 ) return 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 2; ++j ) {
			if ( i == j ) { if ( e3[i + 2 * j] != 1 || e3({i, j}) != 1 ) return 1 ; }
			else { if ( e3[i + 2 * j] != 0 || e3({i, j}) != 0 ) return 1 ; }
    	}
    }

	std::cout << "eye construction" << '\n' ;
	std::cout << "glas3::eye<double, 1> e4({3}) ;" << '\n' ;
	glas3::eye<double, 1> e4({3}) ;
	std::cout << e4 << '\n' << '\n' ;
	if ( e4.size() != 3 || e4.ndof() != 0 || e4.shape().size() != 1 || e4.shape()[0] != 3 || e4.shape().shape().size() != 1 ) return 1 ;
    for ( i = 0; i < 3; ++i ) {
    	if ( i == 0 ) { if ( e4[i] != 1 || e4({i}) != 1 ) return 1 ; }
    	else { if ( e4[i] != 0 || e4({i}) != 0 ) return 1 ; }
    }

	std::cout << "eye construction" << '\n' ;
	std::cout << "glas3::eye<double, 0> e5({}) ;" << '\n' ;
	glas3::eye<double, 0> e5({}) ;
	std::cout << e5 << '\n' << '\n' ;
	if ( e5.size() != 1 || e5.ndof() != 0 || e5.shape().size() != 0 || e5.shape().shape().size() != 1 ) return 1 ;
	if ( e5[0] != 1 ) return 1 ;

	std::cout << "eye construction" << '\n' ;
	std::cout << "glas3::eye<double, 0> e6(5, 0) ;" << '\n' ;
	glas3::eye<double, 0> e6(5, 0) ;
	std::cout << e6 << '\n' << '\n' ;
	if ( e6.size() != 1 || e6.ndof() != 0 || e6.shape().size() != 0 || e6.shape().shape().size() != 1 ) return 1 ;
	if ( e6[0] != 1 ) return 1 ;

	std::cout << "eye construction" << '\n' ;
	std::cout << "glas3::eye<double, 0> e7 ;" << '\n' ;
	glas3::eye<double, 0> e7 ;
	std::cout << e7 << '\n' << '\n' ;
	if ( e7.size() != 1 || e7.ndof() != 0 || e7.shape().size() != 0 || e7.shape().shape().size() != 1 ) return 1 ;
	if ( e7[0] != 1 ) return 1 ;

	std::cout << '\n' ;
	return 0;
}
