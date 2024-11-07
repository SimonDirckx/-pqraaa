#include <glas3/array/type/empty_array.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/linear_index_selection.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j ;

	std::cout << "glas3::dense_vector<int> v1 = {1, 2, 3} ;" << '\n' ;
	std::cout << "glas3::dense_array<int> a1({0, 1, 0, 2, 1, 0}, {2, 3}) ;" << '\n' ;
	std::cout << "auto v1_a1 = v1[a1] ;" << '\n' ;
	glas3::dense_vector<int> v1 = {1, 2, 3} ;
	glas3::dense_array<int> a1({0, 1, 0, 2, 1, 0}, {2, 3}) ;
 	auto v1_a1 = v1[a1] ;

	std::cout << "linear index array selection of indirect_array using Array (shallow copy)" << '\n' ;
    std::cout << "glas3::dense_array<int> a2({0, 1, 5, 2, 3, 0}, {3, 2}) ;" << '\n' ;
    std::cout << "auto v1_a1_a2 = v1_a1[a2] ;" << '\n' ;
    glas3::dense_array<int> a2({0, 1, 5, 2, 3, 0}, {3, 2}) ;
    auto v1_a1_a2 = v1_a1[a2] ;
    if ( v1_a1_a2.size() != 6 || v1_a1_a2.ndof() != v1_a1.ndof() + a2.ndof() || v1_a1_a2.shape().size() != 2 || v1_a1_a2.shape()[0] != 3 || v1_a1_a2.shape()[1] != 2 || v1_a1_a2.shape().shape().size() != 1 ) return 1 ;
    v1[0] -= 1 ;
    a1[1] += 1 ;
    a2[2] -= 1 ;
    for ( i = 0; i < 3; ++i ) {
    	for ( j = 0; j < 2; ++j ) {
    		if ( v1_a1_a2[i + 3 * j] != v1_a1[a2[i + 3 * j]] || v1_a1_a2({i, j}) != v1_a1[a2({i, j})] ) return 1 ;//
    	}
    }

    std::cout << "linear index array selection of indirect_array using std::initializer_list (shallow copy)" << '\n' ;
    std::cout << "auto v1_a1_l = v1_a1[{1, 2, 0, 1, 5}] ;" << '\n' ;
    glas3::dense_vector<int> v2 = {1, 2, 0, 1, 5}  ;
    auto v1_a1_l = v1_a1[{1, 2, 0, 1, 5}] ;
    if ( v1_a1_l.size() != 5 || v1_a1_l.ndof() != v1_a1.ndof() + 5 || v1_a1_l.shape().size() != 1 || v1_a1_l.shape()[0] != 5 || v1_a1_l.shape().shape().size() != 0 ) return 1 ;
    v1[0] -= 1 ;
    for ( i = 0; i < 5; ++i ) {
    	if ( v1_a1_l[i] != v1_a1[v2[i]] || v1_a1_l({i}) != v1_a1[v2({i})] ) return 1 ;
    }

	std::cout << "linear index array selection of indirect_array using empty std::initializer_list (shallow copy)" << '\n' ;
	std::cout << "auto v1_a1_el = v1_a1[{}] ;" << '\n' ;
	auto v1_a1_el = v1_a1[glas3::empty_array()] ;
    if ( v1_a1_el.size() != 0 || v1_a1_el.ndof() != v1_a1.ndof() || v1_a1_el.shape().size() != 0 || v1_a1_el.shape().shape().size() != 0 ) return 1 ;

    std::cout << '\n' ;

	return 0;

}
