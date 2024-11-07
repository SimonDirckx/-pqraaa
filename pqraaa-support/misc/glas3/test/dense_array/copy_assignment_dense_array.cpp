#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <iostream>
#include <vector>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j ;

	std::cout << "copy_assignment of empty dense_array from dense_array, std::initializer_list, Container, Array (deep copy)" << '\n' ;
	std::cout << "glas3::dense_array<double> a1 ;" << '\n' ;
	std::cout << "glas3::dense_array<double> a2 ;" << '\n' ;
	std::cout << "glas3::dense_scalar<double> s1 ;" << '\n' ;
	std::cout << "a1 = a2 ;" << '\n' ;
	std::cout << "a1 = {} ;" << '\n' ;
	std::cout << "a1 = s1 ;" << '\n' ;
	std::cout << "a1 = std::vector<double>() ;" << '\n' ;
	glas3::dense_array<double> a1 ;
	glas3::dense_array<double> a2 ;
	glas3::dense_scalar<double> s1 ;

	a1 = a2 ;
	if ( a1.size() != 0 || a1.ndof() != 0 || a1.shape().size() != 0 || a1.shape().shape().size() != 1 ) return 1 ;

	a1 = {} ;
	if ( a1.size() != 0 || a1.ndof() != 0 || a1.shape().size() != 0 || a1.shape().shape().size() != 1 ) return 1 ;

	a1 = s1 ;
	if ( a1.size() != 0 || a1.ndof() != 0 || a1.shape().size() != 0 || a1.shape().shape().size() != 1 ) return 1 ;

	a1 = std::vector<double>() ;
	if ( a1.size() != 0 || a1.ndof() != 0 || a1.shape().size() != 0 || a1.shape().shape().size() != 1 ) return 1 ;

	std::cout << "copy_assignment of dense_array from dense_array, value, std::initializer_list, Container, Array (deep copy)" << '\n' ;
	std::cout << "glas3::dense_array<double> a3 = 1 ;" << '\n' ;
	std::cout << "glas3::dense_array<double> a4({1, 2, 3, 4, 5, 6}, {2, 3}) ;" << '\n' ;
	std::cout << "glas3::dense_array<double> a5({7, 8, 9, 10, 11, 12}, {3, 2}) ;" << '\n' ;
	std::cout << "glas3::dense_vector<double> v1 = {13, 14, 15, 16, 17, 18} ;" << '\n' ;
	std::cout << "a3 = 2 ;" << '\n' ;
	std::cout << "a4 = a5 ;" << '\n' ;
	std::cout << "a4 = {19, 20, 21, 22, 23, 24} ;" << '\n' ;
	std::cout << "a4 = v1 ;" << '\n' ;
	std::cout << "a4 = std::vector<double>({19, 20, 21, 22, 23, 24}) ;" << '\n' ;
	glas3::dense_array<double> a3 = 1 ;
	glas3::dense_array<double> a4({1, 2, 3, 4, 5, 6}, {2, 3}) ;
	glas3::dense_array<double> a5({7, 8, 9, 10, 11, 12}, {3, 2}) ;
	glas3::dense_vector<double> v1 = {13, 14, 15, 16, 17, 18} ;

	a3 = 2 ;
	if ( a3.size() != 1 || a3.ndof() != 1 || a3.shape().size() != 0 || a3.shape().shape().size() != 1 ) return 1 ;
	if ( a3[0] != 2 ) return 1 ;

	a4 = a5 ;
	if ( a4.size() != 6 || a4.ndof() != 6 || a4.shape().size() != 2 || a4.shape()[0] != 2 || a4.shape()[1] != 3 || a4.shape().shape().size() != 1 ) return 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		auto dum = a4({i, j}) ;
    		a5[i + 2 * j] += 1 ;
    		if ( a4({i, j}) != dum ) return 1 ;
    	}
    }

    glas3::dense_vector<double> v2 = {19, 20, 21, 22, 23, 24} ;
	a4 = {19, 20, 21, 22, 23, 24} ;
	if ( a4.size() != 6 || a4.ndof() != 6 || a4.shape().size() != 2 || a4.shape()[0] != 2 || a4.shape()[1] != 3 || a4.shape().shape().size() != 1 ) return 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		if ( a4({i, j}) != v2[i + 2 * j] ) return 1 ;
    	}
    }

	a4 = v1 ;
	if ( a4.size() != 6 || a4.ndof() != 6 || a4.shape().size() != 2 || a4.shape()[0] != 2 || a4.shape()[1] != 3 || a4.shape().shape().size() != 1 ) return 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		auto dum = a4({i, j}) ;
    		v1[i + 2 * j] += 1 ;
    	    if ( a4({i, j}) != dum ) return 1 ;
    	}
    }

	a4 = v2 ;
	if ( a4.size() != 6 || a4.ndof() != 6 || a4.shape().size() != 2 || a4.shape()[0] != 2 || a4.shape()[1] != 3 || a4.shape().shape().size() != 1 ) return 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		if ( a4({i, j}) != v2[i + 2 * j] ) return 1 ;
    	}
    }

	std::cout << '\n' ;
	return 0;

}
