#include <glas3/array/type/empty_array.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>
#include <glas3/array/dense_array/type/range.hpp>

#include <glas3/array/dense_array/algorithm/linear_index_selection.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j ;

	std::cout << "linear index array selection of range using Array (shallow copy)" << '\n' ;
	std::cout << "glas3::range<> r1(0, 5) ;" << '\n' ;
	std::cout << "glas3::dense_array<int> a1({0, 1, 2, 4, 1, 0}, {2, 3}) ;" << '\n' ;
	std::cout << "auto r1_a1 = r1[a1] ;" << '\n' ;
	glas3::range<> r1(0, 5) ;
	glas3::dense_array<int> a1({0, 3, 2, 4, 1, 0}, {2, 3}) ;
	auto r1_a1 = r1[a1] ;
    if ( r1_a1.size() != 6 || r1_a1.ndof() != r1.ndof() + a1.ndof() || r1_a1.shape().size() != 2 || r1_a1.shape()[0] != 2 || r1_a1.shape()[1] != 3 || r1_a1.shape().shape().size() != 1 ) return 1 ;
    a1[0] += 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		if ( r1_a1[i + 2 * j] != r1[a1[i + 2 * j]] || r1_a1({i, j}) != r1[a1({i, j})] ) return 1 ;
     	}
    }

	std::cout << "linear index array selection of range using std::initializer_list" << '\n' ;
	std::cout << "glas3::range<> r2(0, 5) ;" << '\n' ;
	std::cout << "auto r2_l = r2[{0, 1, 0, 2, 4, 4, 3}] ;" << '\n' ;
	glas3::range<> r2(0, 5) ;
	glas3::dense_vector<int> v1 = {0, 1, 0, 2, 4, 4, 3} ;
	auto r2_l = r2[{0, 1, 0, 2, 4, 4, 3}] ;
    if ( r2_l.size() != 7 || r2_l.ndof() != r2.ndof() + v1.ndof() || r2_l.shape().size() != 1 || r2_l.shape()[0] != 7 || r2_l.shape().shape().size() != 0 ) return 1 ;
    for ( i = 0; i < 4; ++i ) {
    	if ( r2_l[i] != r2[v1[i]] || r2_l({i}) != r2[v1({i})] ) return 1 ;
    }

	std::cout << "linear index array selection of dense_vector using empty std::initializer_list" << '\n' ;
	std::cout << "glas3::range<> r3(0, 5) ;" << '\n' ;
	std::cout << "auto r3_l = r3[{}] ;" << '\n' ;
	glas3::range<> r3(0, 5) ;
	auto r3_l = r3[glas3::empty_array()] ;
    if ( r3_l.size() != 0 || r3_l.ndof() != r3.ndof() || r3_l.shape().size() != 0 || r3_l.shape().shape().size() != 0 ) return 1 ;

    std::cout << '\n' ;

	return 0;

}
