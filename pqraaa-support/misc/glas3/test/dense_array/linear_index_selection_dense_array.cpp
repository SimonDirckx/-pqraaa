#include <glas3/array/type/empty_array.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/linear_index_selection.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i ;

	std::cout << "linear index array selection of dense_array using Array (shallow copy)" << '\n' ;
	std::cout << "glas3::dense_array<int> a1({1, 2, 3, 4, 5, 6}, {2, 3}) ;" << '\n' ;
	std::cout << "glas3::dense_vector<int> v1 = {1, 2, 0, 1, 5} ;" << '\n' ;
	std::cout << "auto a1_v1 = a1[v1] ;" << '\n' ;
	glas3::dense_array<int> a1({1, 2, 3, 4, 5, 6}, {2, 3}) ;
	glas3::dense_vector<int> v1 = {1, 2, 0, 1, 5} ;
 	auto a1_v1 = a1[v1] ;
    if ( a1_v1.size() != 5 || a1_v1.ndof() != 6 + 5 || a1_v1.shape().size() != 1 || a1_v1.shape()[0] != 5 || a1_v1.shape().shape().size() != 0 ) return 1 ;
    a1[0] -= 1 ;
    v1[1] -= 1 ;
    for ( i = 0; i < 5; ++i ) {
    	if ( a1_v1[i] != a1[v1[i]] || a1_v1({i}) != a1[v1({i})] ) return 1 ;
    }

	std::cout << "linear index array selection of dense_array using std::initializer_list (shallow copy)" << '\n' ;
	std::cout << "glas3::dense_array<int> a2({1, 2, 3, 4, 5, 6}, {2, 3}) ;" << '\n' ;
	std::cout << "auto a2_l = a2[{1, 2, 0, 1, 5}] ;" << '\n' ;
	glas3::dense_array<int> a2({1, 2, 3, 4, 5, 6}, {2, 3}) ;
	glas3::dense_vector<int> v2 = {1, 2, 0, 1, 5} ;
	auto a2_l = a2[{1, 2, 0, 1, 5}] ;
    if ( a2_l.size() != 5 || a2_l.ndof() != 6 + 5 || a2_l.shape().size() != 1 || a2_l.shape()[0] != 5 || a2_l.shape().shape().size() != 0 ) return 1 ;
    a2[0] -= 1 ;
    for ( i = 0; i < 5; ++i ) {
    	if ( a2_l[i] != a2[v2[i]] || a2_l({i}) != a2[v2({i})] ) return 1 ;
    }

	std::cout << "linear index array selection of dense_array using empty std::initializer_list (shallow copy)" << '\n' ;
	std::cout << "glas3::dense_array<int> a3({1, 2, 3, 4, 5, 6}, {2, 3}) ;" << '\n' ;
	std::cout << "auto a3_l = a3[{}] ;" << '\n' ;
	glas3::dense_array<int> a3({1, 2, 3, 4, 5, 6}, {2, 3}) ;
	auto a3_l = a3[glas3::empty_array()] ;
    if ( a3_l.size() != 0 || a3_l.ndof() != a3.ndof() || a3_l.shape().size() != 0 || a3_l.shape().shape().size() != 0 ) return 1 ;

    std::cout << '\n' ;

	return 0;

}
