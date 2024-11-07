#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/ttt.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j, k ;

	std::cout << "multi-index array selection of 1D entrywise_operation using empty std::initializer_list (shallow copy) is not implemented" << '\n' ;
//	std::cout << "auto a1_el = a1({{}}) ;" << '\n' ;
//	glas3::dense_array<int> a3({1, 2, 3, 4, 5, 6}, {2, 3}) ;
//	auto a1_el = a1({{}}) ;
//    if ( a1_el.size() != 0 || a1_el.shape().size() != 0 || a1_el.shape().shape().size() != 0 ) return 1 ;

    std::cout << "glas3::dense_vector<int> v1 = {2, 3} ;" << '\n' ;
    std::cout << "glas3::dense_array<int> a1({5, 7, 11, 13, 17, 19}, {3, 2}) ;" << '\n' ;
    std::cout << "auto ttt_v1a1 = glas3::ttt( v1, a1 ) ;" << '\n' ;
 	glas3::dense_vector<int> v1 = {2, 3} ;
 	glas3::dense_array<int> a1({5, 7, 11, 13, 17, 19}, {3, 2}) ;
    auto b1 = glas3::ttt( v1, a1 ) ;

	std::cout << "multi-index array selection of nD dense_array using Array (shallow copy)" << '\n' ;
	std::cout << "glas3::dense_vector<std::size_t> v2 = {1, 2, 0, 1, 2} ;" << '\n' ;
	std::cout << "auto b1_s = b1({1, v2, {1, 0, 1}}) ;" << '\n' ;
	glas3::dense_vector<int> v2 = {1, 2, 0, 1, 2} ;
	glas3::dense_vector<int> v3 = {1, 0, 1} ;
	auto b1_s = b1({1, v2, {1, 0, 1}}) ;
    if ( b1_s.size() != 15 || b1_s.shape().size() != 3 || b1_s.shape()[0] != 1 || b1_s.shape()[1] != 5 || b1_s.shape()[2] != 3 || b1_s.shape().shape().size() != 1 ) return 1 ;
    a1[0] -= 1 ;
    v2[0] -= 1 ;
    for ( i = 0; i < 1; ++i ) {
    	for ( j = 0; j < 5; ++j ) {
    		for ( k = 0; k < 3; ++k ) {
    			if ( b1_s[i + 1 * j + 1 * 5 * k] != b1({1, v2[j], v3[k]}) || b1_s({i, j, k}) != b1({1, v2[j], v3[k]}) ) return 1 ;
    		}
    	}
    }
    std::cout << '\n' ;

	return 0;

}