#include <glas3/array/type/empty_array.hpp>
#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/linear_index_selection.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j ;

	std::cout << "copy_assignment of empty indirect_array from indirect_array, std::initializer_list, Array (deep copy)" << '\n' ;
	std::cout << "glas3::dense_array<double> a1({1, 2, 3, 4, 5, 6}, {2, 3}) ;" << '\n' ;
	std::cout << "glas3::dense_array<double> a2({7, 8, 9, 10, 11, 12}, {3, 2}) ;" << '\n' ;
	std::cout << "glas3::dense_scalar<double> s1 ;" << '\n' ;
	std::cout << "auto a3 = a1[{}] ;" << '\n' ;
	std::cout << "a3 = a2[{}] ;" << '\n' ;
	std::cout << "a3 = {} ;" << '\n' ;
	std::cout << "a3 = s1 ;" << '\n' ;
	glas3::dense_array<double> a1({1, 2, 3, 4, 5, 6}, {2, 3}) ;
	glas3::dense_array<double> a2({7, 8, 9, 10, 11, 12}, {3, 2}) ;
	glas3::dense_scalar<double> s1 ;
	auto a3 = a1[glas3::empty_array()] ;

	if ( a3.size() != 0 || a3.ndof() != a1.ndof() || a3.shape().size() != 0 || a3.shape().shape().size() != 0 ) return 1 ;

	a3 = a2[glas3::empty_array()] ;
	if ( a3.size() != 0 || a2.ndof() != a3.ndof() || a3.shape().size() != 0 || a3.shape().shape().size() != 0 ) return 1 ;

	a3 = {} ;
	if ( a3.size() != 0 || a2.ndof() != a3.ndof() || a3.shape().size() != 0 || a3.shape().shape().size() != 0 ) return 1 ;

	a3 = s1 ;
	if ( a3.size() != 0 || a2.ndof() != a3.ndof() || a3.shape().size() != 0 || a3.shape().shape().size() != 0 ) return 1 ;

	std::cout << "copy_assignment of indirect_array from indirect_array, value, std::initializer_list, Array (deep copy)" << '\n' ;
	std::cout << "glas3::dense_array<double> a4({2, 3, 0, 5}, {2, 2}) ;" << '\n' ;
	std::cout << "auto a5 = a1[a4] ;" << '\n' ;
	std::cout << "glas3::dense_vector<double> v1 = {5, 0, 2, 1} ;" << '\n' ;
	std::cout << "auto a6 = a2[v1] ;" << '\n' ;
	std::cout << "glas3::dense_vector<double> v2 = 1 ;" << '\n' ;
	std::cout << "auto a7 = a1[v2] ;" << '\n' ;
	std::cout << "glas3::dense_vector<double> v3 = {17, 18, 19, 20} ;" << '\n' ;
	std::cout << "a7 = 0 ;" << '\n' ;
	std::cout << "a5 = a6 ;" << '\n' ;
	std::cout << "a5 = {13, 14, 15, 16} ;" << '\n' ;
	std::cout << "a5 = v3 ;" << '\n' ;
	glas3::dense_array<int> a4({2, 3, 0, 5}, {2, 2}) ;
	auto a5 = a1[a4] ;
	glas3::dense_vector<int> v1 = {5, 0, 2, 1} ;
	auto a6 = a2[v1] ;
	glas3::dense_vector<int> v2 = 1 ;
	auto a7 = a1[v2] ;
	glas3::dense_vector<double> v3 = {17, 18, 19, 20} ;

	if ( a7.size() != v2.size() || a7.ndof() != a1.ndof() + v2.ndof() || a7.shape().size() != v2.shape().size() || a7.shape()[0] != v2.shape()[0] || a7.shape().shape().size() != v2.shape().shape().size() ) return 1 ;
	a7 = 0 ;
	if ( a7.size() != v2.size() || a7.ndof() != a1.ndof() + v2.ndof() || a7.shape().size() != v2.shape().size() || a7.shape()[0] != v2.shape()[0] || a7.shape().shape().size() != v2.shape().shape().size() ) return 1 ;
	if ( a7[0] != 0 ) return 1 ;

	a5 = a6 ;
	if ( a5.size() != a4.size() || a5.ndof() != a1.ndof() + a4.ndof() || a5.shape().size() != a4.shape().size() || a5.shape()[0] != a4.shape()[0] || a5.shape()[1] != a4.shape()[1] || a5.shape().shape().size() != a4.shape().shape().size() ) return 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 2; ++j ) {
    		auto dum = a5({i, j}) ;
    		a6[i + 2 * j] += 1 ;
    		if ( a5({i, j}) != dum ) return 1 ;
    	}
    }

    glas3::dense_vector<double> v4 = {13, 14, 15, 16} ;
    a5 = {13, 14, 15, 16} ;
    if ( a5.size() != a4.size() || a5.ndof() != a1.ndof() + a4.ndof() || a5.shape().size() != a4.shape().size() || a5.shape()[0] != a4.shape()[0] || a5.shape()[1] != a4.shape()[1] || a5.shape().shape().size() != a4.shape().shape().size() ) return 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 2; ++j ) {
    		if ( a5({i, j}) != v4[i + 2 * j] ) return 1 ;
    	}
    }

    a5 = v3 ;
    if ( a5.size() != a4.size() || a5.ndof() != a1.ndof() + a4.ndof() || a5.shape().size() != a4.shape().size() || a5.shape()[0] != a4.shape()[0] || a5.shape()[1] != a4.shape()[1] || a5.shape().shape().size() != a4.shape().shape().size() ) return 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 2; ++j ) {
    		auto dum = a5({i, j}) ;
    		v3[i + 2 * j] += 1 ;
    	    if ( a5({i, j}) != dum ) return 1 ;
    	}
    }

	std::cout << '\n' ;
	return 0;

}
