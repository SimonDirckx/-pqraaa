#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_matrix.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j ;

	std::cout << "multi-index array selection of 1D dense_matrix using empty std::initializer_list (shallow copy) is not implemented" << '\n' ;
//	std::cout << "auto a1_el = a1({{}}) ;" << '\n' ;
//	glas3::dense_matrix<int> a3({1, 2, 3, 4, 5, 6}, {2, 3}) ;
//	auto a1_el = a1({{}}) ;
//    if ( a1_el.size() != 0 || a1_el.shape().size() != 0 || a1_el.shape().shape().size() != 0 ) return 1 ;

 	std::cout << "glas3::dense_matrix<int> m1({1, 2, 3, 4, 5, 6}, {2, 3}) ;" << '\n' ;
 	std::cout << "glas3::dense_vector<int> v1 = {1, 2, 0, 1, 2} ;" << '\n' ;
 	glas3::dense_matrix<int> m1({1, 2, 3, 4, 5, 6}, {2, 3}) ;
	glas3::dense_vector<int> v1 = {1, 2, 0, 1, 2} ;

	std::cout << "multi-index array selection of nD dense_matrix using Array (shallow copy)" << '\n' ;
	std::cout << "auto m1_s = m1({1, v1}) ;" << '\n' ;
	auto m1_s = m1({1, v1}) ;
    if ( m1_s.size() != 5 || m1_s.shape().size() != 2 || m1_s.shape()[0] != 1 || m1_s.shape()[1] != 5 || m1_s.shape().shape().size() != 1 ) return 1 ;
    v1[0] -= 1 ;
    m1[1] -= 1 ;
    for ( i = 0; i < 1; ++i ) {
    	for ( j = 0; j < 5; ++j ) {
   			if ( m1_s[i + 1 * j] != m1({1, v1[j]}) || m1_s({i, j}) != m1({1, v1[j]}) ) return 1 ;
    	}
    }

	std::cout << "multi-index array selection of nD dense_matrix using Array (shallow copy)" << '\n' ;
	std::cout << "auto m1_s2 = m1({1, {1, 0, 1}}) ;" << '\n' ;
	glas3::dense_vector<int> v2 = {1, 0, 1} ;
	auto m1_s2 = m1({1, {1, 0, 1}}) ;
    if ( m1_s2.size() != 3 || m1_s2.shape().size() != 2 || m1_s2.shape()[0] != 1 || m1_s2.shape()[1] != 3 || m1_s2.shape().shape().size() != 1 ) return 1 ;
    m1[1] -= 1 ;
    for ( i = 0; i < 1; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
   			if ( m1_s2[i + 1 * j] != m1({1, v2[j]}) || m1_s2({i, j}) != m1({1, v2[j]}) ) return 1 ;
    	}
    }

	std::cout << "multi-index array selection of nD dense_matrix using Array (shallow copy)" << '\n' ;
	std::cout << "auto m1_s3 = m1({{1, 0, 1}, v1}) ;" << '\n' ;
	auto m1_s3 = m1({{1, 0, 1}, v1}) ;
    if ( m1_s3.size() != 15 || m1_s3.shape().size() != 2 || m1_s3.shape()[0] != 3 || m1_s3.shape()[1] != 5 || m1_s3.shape().shape().size() != 1 ) return 1 ;
    v1[0] += 1 ;
    m1[0] -= 1 ;
    for ( i = 0; i < 3; ++i ) {
    	for ( j = 0; j < 5; ++j ) {
    		if ( m1_s3[i + 3 * j] != m1({v2[i], v1[j]}) || m1_s3({i, j}) != m1({v2[i], v1[j]}) ) return 1 ;
    	}
    }

    std::cout << '\n' ;

	return 0;

}
