#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/permute.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j, k ;

	std::cout << "permute operation (shallow copy)" << '\n' ;
	std::cout << "auto permute_a = glas3::permute( glas3::dense_array<int>({1, 2, 3, 4, 5, 6}, {2, 3, 1}), {2, 0, 1} ) ;" << '\n' ;
 	glas3::dense_array<int> a({1, 2, 3, 4, 5, 6}, {2, 3, 1}) ;
 	auto permute_a = glas3::permute( a, {2, 0, 1} ) ;
    if ( permute_a.size() != 6 || permute_a.shape().size() != 3 || permute_a.shape()[0] != 1 || permute_a.shape()[1] != 2 || permute_a.shape()[2] != 3 || permute_a.shape().shape().size() != 1 ) return 1 ;
    a[0] -= 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		for ( k = 0; k < 1; ++k ) {
    			if ( permute_a[k + 1 * i + 1 * 2 * j] != a[i + 2 * j + 2 * 3 * k] || permute_a({k, i, j}) != a({i, j, k}) ) return 1 ;
    		}
    	}
    }

    std::cout << '\n' ;

	return 0;

}
