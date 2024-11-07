#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <iostream>
#include <vector>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j ;

	std::cout << "copy_assignment of empty dense_vector from dense_vector, std::initializer_list, Container, Array (deep copy)" << '\n' ;
	std::cout << "glas3::dense_vector<double> v1 ;" << '\n' ;
	std::cout << "glas3::dense_vector<double> v2 ;" << '\n' ;
	std::cout << "glas3::dense_scalar<double> s1 ;" << '\n' ;
	std::cout << "v1 = v2 ;" << '\n' ;
	std::cout << "v1 = {} ;" << '\n' ;
	std::cout << "v1 = s1 ;" << '\n' ;
	std::cout << "v1 = std::vector<double>() ;" << '\n' ;
	glas3::dense_vector<double> v1 ;
	glas3::dense_vector<double> v2 ;
	glas3::dense_scalar<double> s1 ;

	v1 = v2 ;
	if ( v1.size() != 0 || v1.ndof() != 0 || v1.shape().size() != 1 || v1.shape()[0] != 0 || v1.shape().shape().size() != 0 ) return 1 ;

	v1 = {} ;
	if ( v1.size() != 0 || v1.ndof() != 0 || v1.shape().size() != 1 || v1.shape()[0] != 0 || v1.shape().shape().size() != 0 ) return 1 ;

	v1 = s1 ;
	if ( v1.size() != 0 || v1.ndof() != 0 || v1.shape().size() != 1 || v1.shape()[0] != 0 || v1.shape().shape().size() != 0 ) return 1 ;

	v1 = std::vector<double>() ;
	if ( v1.size() != 0 || v1.ndof() != 0 || v1.shape().size() != 1 || v1.shape()[0] != 0 || v1.shape().shape().size() != 0 ) return 1 ;

	std::cout << "copy_assignment of dense_vector from dense_vector, value, std::initializer_list, Container, Array (deep copy)" << '\n' ;
	std::cout << "glas3::dense_vector<double> v3 = 1 ;" << '\n' ;
	std::cout << "glas3::dense_vector<double> v4 = {1, 2, 3, 4, 5, 6} ;" << '\n' ;
	std::cout << "glas3::dense_vector<double> v5 = {7, 8, 9, 10, 11, 12} ;" << '\n' ;
	std::cout << "glas3::dense_array<double> a1({13, 14, 15, 16, 17, 18}, {2, 3}) ;" << '\n' ;
	std::cout << "v3 = 2 ;" << '\n' ;
	std::cout << "v4 = v5 ;" << '\n' ;
	std::cout << "v4 = {19, 20, 21, 22, 23, 24} ;" << '\n' ;
	std::cout << "v4 = a1 ;" << '\n' ;
	std::cout << "v4 = std::vector<double>({19, 20, 21, 22, 23, 24}) ;" << '\n' ;
	glas3::dense_vector<double> v3 = 1 ;
	glas3::dense_vector<double> v4 = {1, 2, 3, 4, 5, 6} ;
	glas3::dense_vector<double> v5 = {7, 8, 9, 10, 11, 12} ;
	glas3::dense_array<double> a1({13, 14, 15, 16, 17, 18}, {2, 3}) ;

	v3 = 2 ;
	if ( v3.size() != 1 || v3.ndof() != 1 || v3.shape().size() != 1 || v3.shape()[0] != 1 || v3.shape().shape().size() != 0 ) return 1 ;
	if ( v3[0] != 2 ) return 1 ;

	v4 = v5 ;
	if ( v4.size() != 6 || v4.ndof() != 6 || v4.shape().size() != 1 || v4.shape()[0] != 6 || v4.shape().shape().size() != 0 ) return 1 ;
    for ( i = 0; i < 6; ++i ) {
    	auto dum = v4[i] ;
    	v5[i] += 1 ;
    	if ( v4[i] != dum ) return 1 ;
    }

    glas3::dense_vector<double> v6 = {19, 20, 21, 22, 23, 24} ;
	v4 = {19, 20, 21, 22, 23, 24} ;
	if ( v4.size() != 6 || v4.ndof() != 6 || v4.shape().size() != 1 || v4.shape()[0] != 6 || v4.shape().shape().size() != 0 ) return 1 ;
    for ( i = 0; i < 6; ++i ) {
    	if ( v4[i] != v6[i] ) return 1 ;
    }

	v4 = a1 ;
	if ( v4.size() != 6 || v4.ndof() != 6 || v4.shape().size() != 1 || v4.shape()[0] != 6 || v4.shape().shape().size() != 0 ) return 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		auto dum = v4[i + 2 * j] ;
    	    a1({i, j}) += 1 ;
    	    if ( v4[i + 2 * j] != dum ) return 1 ;
    	}
    }

	v4 = v6 ;
	if ( v4.size() != 6 || v4.ndof() != 6 || v4.shape().size() != 1 || v4.shape()[0] != 6 || v4.shape().shape().size() != 0 ) return 1 ;
    for ( i = 0; i < 6; ++i ) {
    	if ( v4[i] != v6[i] ) return 1 ;
    }

	std::cout << '\n' ;
	return 0;

}
