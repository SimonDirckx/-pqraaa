#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j, k ;

	std::cout << "multi-index array selection of 1D dense_array using empty std::initializer_list (shallow copy) is not implemented" << '\n' ;
//	std::cout << "auto a1_el = a1({{}}) ;" << '\n' ;
//	glas3::dense_array<int> a3({1, 2, 3, 4, 5, 6}, {2, 3}) ;
//	auto a1_el = a1({{}}) ;
//    if ( a1_el.size() != 0 || a1_el.shape().size() != 0 || a1_el.shape().shape().size() != 0 ) return 1 ;

 	std::cout << "glas3::dense_array<int> a1({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, {2, 3, 2}) ;" << '\n' ;
 	std::cout << "glas3::dense_vector<int> v1 = {1, 2, 0, 1, 2} ;" << '\n' ;
 	glas3::dense_array<int> a1({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, {2, 3, 2}) ;
	glas3::dense_vector<int> v1({1, 2, 0, 1, 2}) ;

	std::cout << "multi-index array selection of nD dense_array using Array (shallow copy)" << '\n' ;
	std::cout << "auto a1_s = a1({1, v1, {1, 0, 1}}) ;" << '\n' ;
	glas3::dense_vector<int> v2 = {1, 0, 1} ;
	auto a1_s = a1({1, v1, {1, 0, 1}}) ;
    if ( a1_s.size() != 15 || a1_s.shape().size() != 3 || a1_s.shape()[0] != 1 || a1_s.shape()[1] != 5 || a1_s.shape()[2] != 3 || a1_s.shape().shape().size() != 1 ) return 1 ;
    v1[0] += 1 ;
    a1[0] -= 1 ;
    for ( i = 0; i < 1; ++i ) {
    	for ( j = 0; j < 5; ++j ) {
    		for ( k = 0; k < 3; ++k ) {
    			if ( a1_s[i + 1 * j + 1 * 5 * k] != a1({1, v1[j], v2[k]}) || a1_s({i, j, k}) != a1({1, v1[j], v2[k]}) ) return 1 ;
    		}
    	}
    }

    std::cout << '\n' ;

	return 0;

}
