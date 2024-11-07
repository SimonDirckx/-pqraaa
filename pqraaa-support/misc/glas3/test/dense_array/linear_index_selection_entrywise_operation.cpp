#include <glas3/array/type/empty_array.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/linear_index_selection.hpp>

#include <glas3/array/dense_array/algorithm/ops.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j ;

    std::cout << "glas3::dense_array<int> a1({1, 2, 3, 4, 5, 6}, {2, 3}) ;" << '\n' ;
    std::cout << "glas3::dense_vector<int> v1({1, 2, 3, 4, 5, 6}) ;" << '\n' ;
    std::cout << "auto b1 = ( a1 * ( v1 + a1 ) - v1 ) / v1 ;" << '\n' ;
    glas3::dense_array<int> a1({1, 2, 3, 4, 5, 6}, {2, 3}) ;
    glas3::dense_vector<int> v1({1, 2, 3, 4, 5, 6}) ;
    auto b1 = ( a1 * ( v1 + a1 ) - v1 ) / v1 ;

	std::cout << "linear index array selection of entrywise_operation using Array (shallow copy)" << '\n' ;
    std::cout << "glas3::dense_array<int> a2({0, 1, 5, 2, 3, 0}, {3, 2}) ;" << '\n' ;
    std::cout << "auto b1_a2 = b1[a2] ;" << '\n' ;
    glas3::dense_array<int> a2({0, 1, 5, 2, 3, 0}, {3, 2}) ;
    auto b1_a2 = b1[a2] ;
    if ( b1_a2.size() != 6 || b1_a2.ndof() != b1.ndof() + a2.ndof() || b1_a2.shape().size() != 2 || b1_a2.shape()[0] != 3 || b1_a2.shape()[1] != 2 || b1_a2.shape().shape().size() != 1 ) return 1 ;
    a1[0] -= 1 ;
    a2[0] += 1 ;
    for ( i = 0; i < 3; ++i ) {
    	for ( j = 0; j <2; ++j ) {
    		if ( b1_a2[i + 3 * j] != b1[a2[i + 3 * j]] || b1_a2({i, j}) != b1[a2({i, j})] ) return 1 ;
    	}
    }

    std::cout << "linear index array selection of entrywise_operation using std::initializer_list (shallow copy)" << '\n' ;
    std::cout << "auto b1_l = b1[{1, 2, 0, 1, 5}] ;" << '\n' ;
    glas3::dense_vector<int> v2 = {1, 2, 0, 1, 5}  ;
    auto b1_l = b1[{1, 2, 0, 1, 5}] ;
    if ( b1_l.size() != 5 || b1_l.ndof() != b1.ndof() + 5 || b1_l.shape().size() != 1 || b1_l.shape()[0] != 5 || b1_l.shape().shape().size() != 0 ) return 1 ;
    a1[0] -= 1 ;
    for ( i = 0; i < 5; ++i ) {
    	if ( b1_l[i] != b1[v2[i]] || b1_l({i}) != b1[v2({i})] ) return 1 ;
    }

	std::cout << "linear index array selection of entrywise_operation using empty std::initializer_list (shallow copy)" << '\n' ;
	std::cout << "auto b1_el = b1[{}] ;" << '\n' ;
	auto b1_el = b1[glas3::empty_array()] ;
    if ( b1_el.size() != 0 || b1_el.ndof() != b1.ndof() || b1_el.shape().size() != 0 || b1_el.shape().shape().size() != 0 ) return 1 ;

    std::cout << '\n' ;

	return 0;

}
