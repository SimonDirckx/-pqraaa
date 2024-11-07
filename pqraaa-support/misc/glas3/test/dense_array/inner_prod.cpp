#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/inner_prod.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i ;

 	std::cout << "inner_prod operation" << '\n' ;

 	glas3::dense_vector<int> v = {1, 2, 3, 4, 5, 6} ;
 	glas3::dense_array<int> a({-3, 2, -2, 7, -1, 9}, {2, 3, 1}) ;

 	int dum = 0 ;
    for ( i = 0; i < a.size(); ++i ) {
    	dum += a[i] * v[i] ;
    }
    if ( dum != glas3::inner_prod(a, v) ) return 1 ;

    std::cout << '\n' ;

	return 0;

}
