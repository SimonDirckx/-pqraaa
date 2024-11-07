#include <glas3/array/type/empty_array.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_matrix.hpp>

#include <glas3/array/dense_array/algorithm/linear_index_selection.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i ;

	std::cout << "linear index array selection of dense_matrix using Array (shallow copy)" << '\n' ;
	std::cout << "glas3::dense_matrix<int> m1({1, 2, 3, 4, 5, 6}, {2, 3}) ;" << '\n' ;
	std::cout << "glas3::dense_vector<int> v1 = {1, 2, 0, 1, 5} ;" << '\n' ;
	std::cout << "auto m1_v1 = m1[v1] ;" << '\n' ;
	glas3::dense_matrix<int> m1({1, 2, 3, 4, 5, 6}, {2, 3}) ;
	glas3::dense_vector<int> v1 = {1, 2, 0, 1, 5} ;
	auto m1_v1 = m1[v1] ;
    if ( m1_v1.size() != 5 || m1_v1.ndof() != 6 + 5 || m1_v1.shape().size() != 1 || m1_v1.shape()[0] != 5 || m1_v1.shape().shape().size() != 0 ) return 1 ;
    m1[0] -= 1 ;
    v1[0] -= 1 ;
    for ( i = 0; i < 5; ++i ) {
    	if ( m1_v1[i] != m1[v1[i]] || m1_v1({i}) != m1[v1({i})] ) return 1 ;
    }

	std::cout << "linear index array selection of dense_matrix using std::initializer_list (shallow copy)" << '\n' ;
	std::cout << "glas3::dense_matrix<int> m2({1, 2, 3, 4, 5, 6}, {2, 3}) ;" << '\n' ;
	std::cout << "auto m2_l = m2[{1, 2, 0, 1, 5}] ;" << '\n' ;
	glas3::dense_matrix<int> m2({1, 2, 3, 4, 5, 6}, {2, 3}) ;
	glas3::dense_vector<int> v2 = {1, 2, 0, 1, 5} ;
	auto m2_l = m2[{1, 2, 0, 1, 5}] ;
    if ( m2_l.size() != 5 || m2_l.ndof() != 6 + 5 || m2_l.shape().size() != 1 || m2_l.shape()[0] != 5 || m2_l.shape().shape().size() != 0 ) return 1 ;
    m2[0] -= 1 ;
    for ( i = 0; i < 5; ++i ) {
    	if ( m2_l[i] != m2[v2[i]] || m2_l({i}) != m2[v2({i})] ) return 1 ;
    }

	std::cout << "linear index array selection of dense_matrix using empty std::initializer_list (shallow copy)" << '\n' ;
	std::cout << "glas3::dense_matrix<int> m3({1, 2, 3, 4, 5, 6}, {2, 3}) ;" << '\n' ;
	std::cout << "auto m3_l = m3[{}] ;" << '\n' ;
	glas3::dense_matrix<int> m3({1, 2, 3, 4, 5, 6}, {2, 3}) ;
	auto m3_l = m3[glas3::empty_array()] ;
    if ( m3_l.size() != 0 || m3_l.ndof() != m3.ndof() || m3_l.shape().size() != 0 || m3_l.shape().shape().size() != 0 ) return 1 ;

    std::cout << '\n' ;

	return 0;

}
