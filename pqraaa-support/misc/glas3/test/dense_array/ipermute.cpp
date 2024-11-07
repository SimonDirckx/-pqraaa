#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/ipermute.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j, k ;

	std::cout << "ipermute operation (shallow copy)" << '\n' ;
	std::cout << "auto ipermute_a = glas3::ipermute( glas3::dense_array<int>({1, 2, 3, 4, 5, 6}, {2, 3, 1}), {2, 0, 1} ) ;" << '\n' ;
 	glas3::dense_array<int> a({1, 2, 3, 4, 5, 6}, {2, 3, 1}) ;
 	auto ipermute_a = glas3::ipermute( a, {2, 0, 1} ) ;
    if ( ipermute_a.size() != 6 || ipermute_a.shape().size() != 3 || ipermute_a.shape()[0] != 3 || ipermute_a.shape()[1] != 1 || ipermute_a.shape()[2] != 2 || ipermute_a.shape().shape().size() != 1 ) return 1 ;
    a[0] -= 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		for ( k = 0; k < 1; ++k ) {
    			if ( ipermute_a[j + 3 * k + 3 * 1 * i] != a[i + 2 * j + 2 * 3 * k]
    			     || ipermute_a({j, k, i}) != a({i, j, k}) ) return 1 ;
    		}
    	}
    }

	std::cout << "ipermute operation (shallow copy)" << '\n' ;
	std::cout << "auto ipermute_s = glas3::ipermute( glas3::dense_scalar<int>(1), {} ) ;" << '\n' ;
 	glas3::dense_scalar<int> s = 1 ;
 	auto ipermute_s = glas3::ipermute( s, {} ) ;
    if ( ipermute_s.size() != 1 || ipermute_s.shape().size() != 0 || ipermute_s.shape().shape().size() != 1 ) return 1 ;
    s[0] -= 1 ;
    if ( ipermute_s[0] != s[0] ) return 1 ;

    std::cout << '\n' ;

	return 0;

}
