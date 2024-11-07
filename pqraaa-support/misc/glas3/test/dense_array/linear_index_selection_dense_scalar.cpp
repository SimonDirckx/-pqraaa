#include <glas3/array/type/empty_array.hpp>
#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/linear_index_selection.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j ;

	std::cout << "linear index array selection of dense_scalar using Array (shallow copy)" << '\n' ;
	std::cout << "glas3::dense_scalar<int> s1 = 7 ;" << '\n' ;
	std::cout << "glas3::dense_array<int> a1(0, {2, 3}) ;" << '\n' ;
	std::cout << "auto s1_a1 = s1[a1] ;" << '\n' ;
	glas3::dense_scalar<int> s1 = 7 ;
	glas3::dense_array<int> a1(0, {2, 3}) ;
 	auto s1_a1 = s1[a1] ;
    if ( s1_a1.size() != 6 || s1_a1.ndof() != 1 + 6 || s1_a1.shape().size() != 2 || s1_a1.shape()[0] != 2 || s1_a1.shape()[1] != 3 || s1_a1.shape().shape().size() != 1 ) return 1 ;
    s1[0] -= 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		if ( s1_a1[i + 2 * j] != s1[a1[i + 2 * j]] || s1_a1({i, j}) != s1[a1({i, j})] ) return 1 ;
     	}
    }

	std::cout << "linear index array selection of dense_scalar using std::initializer_list (shallow copy)" << '\n' ;
	std::cout << "glas3::dense_scalar<int> s2 = 7 ;" << '\n' ;
	std::cout << "auto s2_l = s2[{0, 0, 0}] ;" << '\n' ;
	glas3::dense_scalar<int> s2 = 7 ;
	glas3::dense_vector<int> v1 = {0, 0, 0} ;
	auto s2_l = s2[{0, 0, 0}] ;
    if ( s2_l.size() != 3 || s2_l.ndof() != 1 + 3 || s2_l.shape().size() != 1 || s2_l.shape()[0] != 3 || s2_l.shape().shape().size() != 0 ) return 1 ;
    s2[0] -= 1 ;
    for ( i = 0; i < 3; ++i ) {
    	if ( s2_l[i] != s2[v1[i]] || s2_l({i}) != s2[v1({i})] ) return 1 ;
    }

	std::cout << "linear index array selection of dense_scalar using empty std::initializer_list (shallow copy)" << '\n' ;
	std::cout << "glas3::dense_scalar<int> s3 = 7 ;" << '\n' ;
	std::cout << "auto s3_l = s3[{}] ;" << '\n' ;
	glas3::dense_scalar<int> s3 = 7 ;
	auto s3_l = s3[glas3::empty_array()] ;
    if ( s3_l.size() != 0 || s3_l.ndof() != s3.ndof() || s3_l.shape().size() != 0 || s3_l.shape().shape().size() != 0 ) return 1 ;

    std::cout << '\n' ;

	return 0;

}
