#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

//#include <glas3/array/dense_array/algorithm/linear_index_selection.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j ;

	std::cout << "multi-index array selection of dense_vector using Array (shallow copy)" << '\n' ;
	std::cout << "glas3::dense_vector<int> v1 = {1, 2, 3} ;" << '\n' ;
	std::cout << "glas3::dense_array<int> a1({0, 1, 0, 2, 1, 0}, {2, 3}) ;" << '\n' ;
	std::cout << "auto v1_a1 = v1({a1}) ;" << '\n' << '\n' ;
	glas3::dense_vector<int> v1 = {1, 2, 3} ;
	glas3::dense_array<int> a1({0, 1, 0, 2, 1, 0}, {2, 3}) ;
	auto v1_a1 = v1({a1}) ;
	if ( v1_a1.size() != 6 || v1_a1.ndof() != v1.ndof() + a1.ndof() || v1_a1.shape().size() != 1 || v1_a1.shape()[0] != 6 || v1_a1.shape().shape().size() != 0 ) return 1 ;
    a1[1] += 1 ;
	v1[0] -= 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		if ( v1_a1[i + 2 * j] != v1[a1[i + 2 * j]] || v1_a1({i + 2 * j}) != v1[a1({i, j})] ) return 1 ;
     	}
    }

	std::cout << "multi-index array selection of dense_vector using std::initializer_list (shallow copy)" << '\n' ;
	std::cout << "glas3::dense_vector<int> v2 = {1, 2, 3} ;" << '\n' ;
	std::cout << "auto v2_l = v2({{0, 1, 0, 2}}) ;" << '\n' << '\n' ;
	glas3::dense_vector<int> v2 = {1, 2, 3} ;
	glas3::dense_vector<int> v3 = {0, 1, 0, 2} ;
	auto v2_l = v2({{0, 1, 0, 2}}) ;
    if ( v2_l.size() != 4 || v2_l.ndof() != v2.ndof() + v3.ndof() || v2_l.shape().size() != 1 || v2_l.shape()[0] != 4 || v2_l.shape().shape().size() != 0 ) return 1 ;
    v2[0] -= 1 ;
    for ( i = 0; i < 4; ++i ) {
    	if ( v2_l[i] != v2[v3[i]] || v2_l({i}) != v2[v3({i})] ) return 1 ;
    }

	std::cout << "multi-index array selection of dense_vector using empty std::initializer_list (shallow copy) is not implemented" << '\n' ;
//	std::cout << "glas3::dense_vector<int> v4 = {1, 2, 3} ;" << '\n' ;
//	std::cout << "auto v4_l = v4({{}}) ;" << '\n' ;
//	glas3::dense_vector<int> v4 = {1, 2, 3} ;
//	auto v4_l = v4({{}}) ;
//    if ( v4_l.size() != 0 || v4_l.shape().size() != 0 || v4_l.shape().shape().size() != 0 ) return 1 ;

    std::cout << '\n' ;

	return 0;

}
