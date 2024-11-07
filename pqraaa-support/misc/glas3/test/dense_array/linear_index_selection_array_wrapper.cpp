#include <glas3/array/type/empty_array.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>
#include <glas3/array/type/array_wrapper.hpp>

#include <glas3/array/dense_array/algorithm/linear_index_selection.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i ;

	std::cout << "glas3::array_wrapper<int> w({1, 2, 3, 4, 5, 6}, {2, 3}) ;" << '\n' << '\n' ;
	glas3::dense_array<double> a( {1, 2, 3, 4, 5, 6}, {2, 3} ) ;
 	glas3::array_wrapper<int> w( a ) ;

 	std::cout << "linear index array selection of array_wrapper using Array (shallow copy)" << '\n' ;
	std::cout << "glas3::dense_vector<int> v1 = {1, 2, 0, 1, 5} ;" << '\n' ;
	std::cout << "auto w_v1 = w[v1] ;" << '\n' ;
	glas3::dense_vector<int> v1 = {1, 2, 0, 1, 5} ;
 	auto w_v1 = w[v1] ;
    if ( w_v1.size() != 5 || w_v1.ndof() != 6 + 5 || w_v1.shape().size() != 1 || w_v1.shape()[0] != 5 || w_v1.shape().shape().size() != 0 ) return 1 ;
    a[0] -= 1 ;
    v1[1] -= 1 ;
    for ( i = 0; i < 5; ++i ) {
    	if ( w_v1[i] != w[v1[i]] || w_v1({i}) != w[v1({i})] ) return 1 ;
    }

	std::cout << "linear index array selection of array_wrapper using std::initializer_list (shallow copy)" << '\n' ;
	std::cout << "auto w_l = w[{1, 2, 0, 1, 5}] ;" << '\n' ;
	glas3::dense_vector<int> v2 = {1, 2, 0, 1, 5} ;
	auto w_l = w[{1, 2, 0, 1, 5}] ;
    if ( w_l.size() != 5 || w_l.ndof() != 6 + 5 || w_l.shape().size() != 1 || w_l.shape()[0] != 5 || w_l.shape().shape().size() != 0 ) return 1 ;
    a[2] -= 1 ;
    for ( i = 0; i < 5; ++i ) {
    	if ( w_l[i] != w[v2[i]] || w_l({i}) != w[v2({i})] ) return 1 ;
    }

	std::cout << "linear index array selection of array_wrapper using empty std::initializer_list (shallow copy)" << '\n' ;
	std::cout << "auto w_e = w[{}] ;" << '\n' ;
	auto w_e = w[glas3::empty_array()] ;
    if ( w_e.size() != 0 || w_e.ndof() != w.ndof() || w_e.shape().size() != 0 || w_e.shape().shape().size() != 0 ) return 1 ;

    std::cout << '\n' ;

	return 0;

}
