#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_matrix.hpp>

#include <iostream>
#include <vector>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j ;

	std::cout << "copy_assignment of empty dense_matrix from dense_matrix, std::initializer_list, Container, Array (deep copy)" << '\n' ;
	std::cout << "glas3::dense_matrix<double> m1 ;" << '\n' ;
	std::cout << "glas3::dense_matrix<double> m2 ;" << '\n' ;
	std::cout << "glas3::dense_scalar<double> s1 ;" << '\n' ;
	std::cout << "m1 = m2 ;" << '\n' ;
	std::cout << "m1 = {} ;" << '\n' ;
	std::cout << "m1 = s1 ;" << '\n' ;
	std::cout << "m1 = std::vector<double>() ;" << '\n' ;
	glas3::dense_matrix<double> m1 ;
	glas3::dense_matrix<double> m2 ;
	glas3::dense_scalar<double> s1 ;

	m1 = m2 ;
	if ( m1.size() != 0 || m1.ndof() != 0 || m1.shape().size() != 2 || m1.shape()[0] != 0 || m1.shape()[1] != 0 || m1.shape().shape().size() != 1 ) return 1 ;

	m1 = {} ;
	if ( m1.size() != 0 || m1.ndof() != 0 || m1.shape().size() != 2 || m1.shape()[0] != 0 || m1.shape()[1] != 0 || m1.shape().shape().size() != 1 ) return 1 ;

	m1 = s1 ;
	if ( m1.size() != 0 || m1.ndof() != 0 || m1.shape().size() != 2 || m1.shape()[0] != 0 || m1.shape()[1] != 0 || m1.shape().shape().size() != 1 ) return 1 ;

	m1 = std::vector<double>() ;
	if ( m1.size() != 0 || m1.ndof() != 0 || m1.shape().size() != 2 || m1.shape()[0] != 0 || m1.shape()[1] != 0 || m1.shape().shape().size() != 1 ) return 1 ;

	std::cout << "copy_assignment of dense_matrix from dense_matrix, value, std::initializer_list, Container, Array (deep copy)" << '\n' ;
	std::cout << "glas3::dense_matrix<double> m3 = 1 ;" << '\n' ;
	std::cout << "glas3::dense_matrix<double> m4({1, 2, 3, 4, 5, 6}, {2, 3}) ;" << '\n' ;
	std::cout << "glas3::dense_matrix<double> m5({7, 8, 9, 10, 11, 12}, {3, 2}) ;" << '\n' ;
	std::cout << "glas3::dense_vector<double> v1 = {13, 14, 15, 16, 17, 18} ;" << '\n' ;
	std::cout << "m3 = 2 ;" << '\n' ;
	std::cout << "m4 = m5 ;" << '\n' ;
	std::cout << "m4 = {19, 20, 21, 22, 23, 24} ;" << '\n' ;
	std::cout << "m4 = v1 ;" << '\n' ;
	std::cout << "m4 = glas3::dense_vector<double>({19, 20, 21, 22, 23, 24}) ;" << '\n' ;
	glas3::dense_matrix<double> m3 = 1 ;
	glas3::dense_matrix<double> m4({1, 2, 3, 4, 5, 6}, {2, 3}) ;
	glas3::dense_matrix<double> m5({7, 8, 9, 10, 11, 12}, {3, 2}) ;
	glas3::dense_vector<double> v1 = {13, 14, 15, 16, 17, 18} ;

	m3 = 2 ;
	if ( m3.size() != 1 || m3.ndof() != 1 || m3.shape().size() != 2 || m3.shape()[0] != 1 || m3.shape()[1] != 1 || m3.shape().shape().size() != 1 ) return 1 ;
	if ( m3[0] != 2 ) return 1 ;

	m4 = m5 ;
	if ( m4.size() != 6 || m4.ndof() != 6 || m4.shape().size() != 2 || m4.shape()[0] != 2 || m4.shape()[1] != 3 || m4.shape().shape().size() != 1 ) return 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		auto dum = m4({i, j}) ;
    		m5[i + 2 * j] += 1 ;
    		if ( m4({i, j}) != dum ) return 1 ;
    	}
    }

    glas3::dense_vector<double> v2 = {19, 20, 21, 22, 23, 24} ;
	m4 = {19, 20, 21, 22, 23, 24} ;
	if ( m4.size() != 6 || m4.ndof() != 6 || m4.shape().size() != 2 || m4.shape()[0] != 2 || m4.shape()[1] != 3 || m4.shape().shape().size() != 1 ) return 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		if ( m4({i, j}) != v2[i + 2 * j] ) return 1 ;
    	}
    }

	m4 = v1 ;
	if ( m4.size() != 6 || m4.ndof() != 6 || m4.shape().size() != 2 || m4.shape()[0] != 2 || m4.shape()[1] != 3 || m4.shape().shape().size() != 1 ) return 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		auto dum = m4({i, j}) ;
    		v1[i + 2 * j] += 1 ;
    	    if ( m4({i, j}) != dum ) return 1 ;
    	}
    }

    m4 = v2 ;
	if ( m4.size() != 6 || m4.ndof() != 6 || m4.shape().size() != 2 || m4.shape()[0] != 2 || m4.shape()[1] != 3 || m4.shape().shape().size() != 1 ) return 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		if ( m4({i, j}) != v2[i + 2 * j] ) return 1 ;
    	}
    }

	std::cout << '\n' ;
	return 0;

}
