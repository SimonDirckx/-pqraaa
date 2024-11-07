#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/ops.hpp>
#include <glas3/array/dense_array/algorithm/any.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::cout << "any operation" << '\n' ;
 	glas3::dense_vector<int> v = {1, 2, 3, 4, 5, 6} ;
 	glas3::dense_array<int> a({1, 2, 3, 4, 5, 7}, {2, 3, 1}) ;

    if ( glas3::any( v != v ) || !glas3::any( v == a ) ) return 1 ;

    std::cout << '\n' ;

	return 0;

}
