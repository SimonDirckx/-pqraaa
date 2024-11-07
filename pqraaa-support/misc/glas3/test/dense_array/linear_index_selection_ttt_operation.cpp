#include <glas3/array/type/empty_array.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/linear_index_selection.hpp>

#include <glas3/array/dense_array/algorithm/ttt.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j ;

    std::cout << "glas3::dense_vector<int> v1 = {2, 3} ;" << '\n' ;
    std::cout << "glas3::dense_array<int> a1({5, 7, 11, 13}, {2, 2}) ;" << '\n' ;
    std::cout << "auto ttt_v1a1 = glas3::ttt( v1, a1 ) ;" << '\n' ;
 	glas3::dense_vector<int> v1 = {2, 3} ;
 	glas3::dense_array<int> a1({5, 7, 11, 13}, {2, 2}) ;
    auto ttt_v1a1 = glas3::ttt( v1, a1 ) ;

	std::cout << "linear index array selection of ttt_operation using Array (shallow copy)" << '\n' ;
    std::cout << "glas3::dense_array<int> a2({0, 1, 5, 2, 3, 0}, {3, 2}) ;" << '\n' ;
    std::cout << "auto ttt_v1a1_a2 = ttt_v1a1[a2] ;" << '\n' ;
    glas3::dense_array<int> a2({0, 1, 5, 2, 3, 0}, {3, 2}) ;
    auto ttt_v1a1_a2 = ttt_v1a1[a2] ;
    if ( ttt_v1a1_a2.size() != 6 || ttt_v1a1_a2.ndof() != ttt_v1a1.ndof() + a2.ndof() || ttt_v1a1_a2.shape().size() != 2 || ttt_v1a1_a2.shape()[0] != 3 || ttt_v1a1_a2.shape()[1] != 2 || ttt_v1a1_a2.shape().shape().size() != 1 ) return 1 ;
    a1[0] -= 1 ;
    a2[0] += 1 ;
    for ( i = 0; i < 3; ++i ) {
    	for ( j = 0; j <2; ++j ) {
    		if ( ttt_v1a1_a2[i + 3 * j] != ttt_v1a1[a2[i + 3 * j]] || ttt_v1a1_a2({i, j}) != ttt_v1a1[a2({i, j})] ) return 1 ;
    	}
    }

    std::cout << "linear index array selection of ttt_operation using std::initializer_list (shallow copy)" << '\n' ;
    std::cout << "auto ttt_v1a1_l = ttt_v1a1[{1, 2, 0, 1, 5}] ;" << '\n' ;
    glas3::dense_vector<int> v2 = {1, 2, 0, 1, 5}  ;
    auto ttt_v1a1_l = ttt_v1a1[{1, 2, 0, 1, 5}] ;
    if ( ttt_v1a1_l.size() != 5 || ttt_v1a1_l.ndof() != ttt_v1a1.ndof() + 5 || ttt_v1a1_l.shape().size() != 1 || ttt_v1a1_l.shape()[0] != 5 || ttt_v1a1_l.shape().shape().size() != 0 ) return 1 ;
    a1[0] -= 1 ;
    for ( i = 0; i < 5; ++i ) {
    	if ( ttt_v1a1_l[i] != ttt_v1a1[v2[i]] || ttt_v1a1_l({i}) != ttt_v1a1[v2({i})] ) return 1 ;
    }

	std::cout << "linear index array selection of ttt_operation using empty std::initializer_list (shallow copy)" << '\n' ;
	std::cout << "auto ttt_v1a1_el = ttt_v1a1[{}] ;" << '\n' ;
	auto ttt_v1a1_el = ttt_v1a1[glas3::empty_array()] ;
    if ( ttt_v1a1_el.size() != 0 || ttt_v1a1_el.ndof() != ttt_v1a1.ndof() || ttt_v1a1_el.shape().size() != 0 || ttt_v1a1_el.shape().shape().size() != 0 ) return 1 ;

    std::cout << '\n' ;

	return 0;

}
