#include <glas3/array/dense_array/container/dense_array.hpp>
#include <glas3/array/type/array_wrapper.hpp>

#include <iostream>
#include <vector>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j ;

	std::cout << "empty array_wrapper construction" << '\n' ;
	std::cout << "glas3::array_wrapper<double> w1 ;" << '\n' ;
	glas3::array_wrapper<double> w1 ;
	if ( w1.size() != 0 || w1.ndof() != 0 || w1.shape().size() != 0 || w1.shape().shape().size() != 1 ) return 1 ;

	std::cout << "array_wrapper construction (deep copy) from value" << '\n' ;
	std::cout << "glas3::array_wrapper<int> w2(7) ;" << '\n' ;
	glas3::array_wrapper<int> w2(7) ;
	if ( w2.size() != 1 || w2.ndof() != 1 || w2.shape().size() != 0 || w2.shape().shape().size() != 1 ) return 1 ;
	if ( w2[0] != 7 ) return 1 ;

	std::cout << "array_wrapper construction (deep copy) from std::initializer_list" << '\n' ;
	std::cout << "glas3::array_wrapper<int> w3({1, 2}) ;" << '\n' ;
	glas3::array_wrapper<int> w3({1, 2}) ;
	if ( w3.size() != 2 || w3.ndof() != 2 || w3.shape().size() != 1 || w3.shape()[0] != 2 || w3.shape().shape().size() != 1 ) return 1 ;
    if ( w3[0] != 1 || w3[1] != 2 || w3({0}) != 1 || w3({1}) != 2 ) return 1 ;

	std::cout << "array_wrapper construction (shallow copy) from array_wrapper" << '\n' ;
	std::cout << "glas3::array_wrapper<double> w4( a ) ;" << '\n' ;
	glas3::dense_array<double> a( {1, 2, 3, 4, 5, 6}, {2, 3} ) ;
	glas3::array_wrapper<double> w4 = ( a ) ;
	if ( w4.size() != a.size() || w4.ndof() != a.ndof() || w4.shape().size() != a.shape().size() || w4.shape()[0] != a.shape()[0] || w4.shape()[1] != a.shape()[1] || w4.shape().shape().size() != 1 ) return 1 ;
	a[0] -= 1 ;
	for ( i = 0; i < 2; ++i ) {
		for ( j = 0; j < 3; ++j ) {
			if ( w4[i + 2 * j] != a[i + 2 * j] || w4({i, j}) != a({i, j}) ) return 1 ;
		}
	}

	std::cout << '\n' ;

	return 0;
}
