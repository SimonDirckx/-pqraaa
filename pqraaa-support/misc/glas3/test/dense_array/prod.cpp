#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/prod.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i ;

 	std::cout << "prod operation" << '\n' ;

 	glas3::dense_array<int> a({-3, 2, -2, 7, -1, 9}, {2, 3, 1}) ;

 	int dum = 1 ;
    for ( i = 0; i < a.size(); ++i ) {
    	dum *= a[i] ;
    }
    if ( dum != glas3::prod(a) ) return 1 ;

    std::cout << '\n' ;

	return 0;

}
