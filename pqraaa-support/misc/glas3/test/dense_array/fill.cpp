#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/fill.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i ;

 	std::cout << "fill operation" << '\n' ;
 	glas3::dense_array<int> a({1, 2, 3, 4, 5, 6}, {2, 3, 1}) ;
    glas3::fill(a, 0) ;

    for ( i = 0; i < a.size(); ++i ) {
    	if ( a[0] != 0 ) return 1 ;
    }

    std::cout << '\n' ;

	return 0;

}
