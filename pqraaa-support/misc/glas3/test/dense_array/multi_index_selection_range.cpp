#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>
#include <glas3/array/dense_array/type/range.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j ;

	std::cout << "multi-index array selection of range using Array" << '\n' ;
	std::cout << "glas3::range<> r1(5, 10) ;" << '\n' ;
	std::cout << "glas3::dense_array<int> a1({0, 3, 0, 2, 1, 4}, {2, 3}) ;" << '\n' ;
	std::cout << "auto r1_a1 = r1({a1}) ;" << '\n' ;
	glas3::range<> r1(5, 10) ;
	glas3::dense_array<int> a1({0, 3, 0, 2, 1, 4}, {2, 3}) ;
 	auto r1_a1 = r1({a1}) ;
    if ( r1_a1.size() != 6 || r1_a1.shape().size() != 1 || r1_a1.shape()[0] != 6 || r1_a1.shape().shape().size() != 0 ) return 1 ;
    a1[0] += 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		std::cout << r1_a1[i + 2 * j] << ' ' ;
    		if ( r1_a1[i + 2 * j] != r1[a1[i + 2 * j]] || r1_a1({i + 2 * j}) != r1({a1({i, j})}) ) return 1 ;
     	}
    }
    std::cout << '\n' ;

	std::cout << "multi-index array selection of range using std::initializer_list" << '\n' ;
	std::cout << "glas3::range<> r2(5, 10) ;" << '\n' ;
	std::cout << "auto r2_l = r2({{4, 1, 0, 2, 4, 3, 2}}) ;" << '\n' ;
	glas3::range<> r2(5, 10) ;
	glas3::dense_vector<int> v1 = {4, 1, 0, 2, 4, 3, 2} ;
	auto r2_l = r2({{4, 1, 0, 2, 4, 3, 2}}) ;
    if ( r2_l.size() != 7 || r2_l.shape().size() != 1 || r2_l.shape()[0] != 7 || r2_l.shape().shape().size() != 0 ) return 1 ;
    for ( i = 0; i < 7; ++i ) {
    	if ( r2_l[i] != r2[v1[i]] || r2_l({i}) != r2({v1({i})}) ) return 1 ;
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
