#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/iostream.hpp>
#include <glas3/array/dense_array/algorithm/transpose.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::cout << "iostream" << '\n' ;
 	glas3::dense_vector<int> v = {1, 2, 3, 4, 5, 6} ;
 	glas3::dense_array<int> a({1, 2, 3, 4, 5, 7}, {2, 3, 1}) ;

    std::cout << v << '\n' ;
    std::cout << a << '\n' ;
    std::cout << glas3::transpose(a) << '\n' ;

    std::cout << '\n' ;

	return 0;

}
