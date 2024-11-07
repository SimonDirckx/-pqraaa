#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/transpose.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j ;

	std::cout << "transpose operation (shallow copy)" << '\n' ;
	std::cout << "auto transpose_a = glas3::transpose( glas3::dense_array<int>({1, 2, 3, 4, 5, 6}, {2, 3}) ) ;" << '\n' ;
 	glas3::dense_array<int> a({1, 2, 3, 4, 5, 6}, {2, 3}) ;
 	auto transpose_a = glas3::transpose( a ) ;
    if ( transpose_a.size() != 6 || transpose_a.shape().size() != 2 || transpose_a.shape()[0] != 3 || transpose_a.shape()[1] != 2 || transpose_a.shape().shape().size() != 1 ) return 1 ;
    a[0] -= 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		if ( transpose_a[j + 3 * i] != a[i + 2 * j] || transpose_a({j, i}) != a({i, j}) ) return 1 ;
    	}
    }

//	std::cout << "transpose operation (shallow copy)" << '\n' ;
//	std::cout << "auto transpose_s = glas3::transpose( glas3::dense_scalar<int>(1) ) ;" << '\n' ;
// 	glas3::dense_scalar<int> s = 1 ;
// 	auto transpose_s = glas3::transpose( s ) ;
//    if ( transpose_s.size() != 1 || transpose_s.shape().size() != 0 || transpose_s.shape().shape().size() != 1 ) return 1 ;
//    s[0] -= 1 ;
//    if ( transpose_s[0] != s[0] ) return 1 ;

    std::cout << '\n' ;

	return 0;

}
