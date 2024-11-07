#include <glas3/array/dense_array/container/dense_array.hpp>

#include <iostream>
#include <vector>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j ;

	std::cout << "empty dense_array construction" << '\n' ;
	std::cout << "glas3::dense_array<double> a1 ;" << '\n' ;
	glas3::dense_array<double> a1 ;
	if ( a1.size() != 0 || a1.ndof() != 0 || a1.shape().size() != 0 || a1.shape().shape().size() != 1 ) return 1 ;

	std::cout << "dense_array construction (deep copy) from value" << '\n' ;
	std::cout << "glas3::dense_array<int> a2(7) ;" << '\n' ;
	glas3::dense_array<int> a2(7) ;
	if ( a2.size() != 1 || a2.ndof() != 1 || a2.shape().size() != 0 || a2.shape().shape().size() != 1 ) return 1 ;
	if ( a2[0] != 7 ) return 1 ;

	std::cout << "dense_array construction (deep copy) from std::initializer_list" << '\n' ;
	std::cout << "glas3::dense_array<int> a3({1, 2}) ;" << '\n' ;
	glas3::dense_array<int> a3({1, 2}) ;
	if ( a3.size() != 2 || a3.ndof() != 2 || a3.shape().size() != 1 || a3.shape()[0] != 2 || a3.shape().shape().size() != 1 ) return 1 ;
    if ( a3[0] != 1 || a3[1] != 2 || a3({0}) != 1 || a3({1}) != 2 ) return 1 ;

	std::cout << "dense_array construction (deep copy) from Container std::vector" << '\n' ;
	std::cout << "glas3::dense_array<int> a8(std::vector<int>({1, 2})) ;" << '\n' ;
	glas3::dense_array<int> a8(std::vector<int>({1, 2})) ;
	if ( a8.size() != 2 || a8.ndof() != 2 || a8.shape().size() != 1 || a8.shape()[0] != 2 || a8.shape().shape().size() != 1 ) return 1 ;
    if ( a8[0] != 1 || a8[1] != 2 || a8({0}) != 1 || a8({1}) != 2 ) return 1 ;

    std::cout << "dense_array construction (deep copy) from value and shape which can be value, std::initializer_list, Container, or Array" << '\n' ;
    std::cout << "glas3::dense_array<int> a4(0, {2, 3}) ;" << '\n' ;
    glas3::dense_array<int> a4(0, {2, 3}) ;
    if ( a4.size() != 6 || a4.ndof() != 6 || a4.shape().size() != 2 || a4.shape()[0] != 2 || a4.shape()[1] != 3 || a4.shape().shape().size() != 1 ) return 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		if ( a4[i + 2 * j] != 0 || a4({i, j}) != 0 ) return 1 ;
    	}
    }

    std::cout << "uninitialized dense_array construction (deep copy) from shape -> boost::make_shared_noinit apparently returns default initialized??" << '\n' ;
    std::cout << "glas3::dense_array<int> a5(glas3::no_init(), {2, 3}) ;" << '\n' ;
    glas3::dense_array<int> a5(glas3::no_init(), {2, 3}) ;
    if ( a5.size() != 6 || a5.ndof() != 6 || a5.shape().size() != 2 || a5.shape()[0] != 2 || a5.shape()[1] != 3 || a5.shape().shape().size() != 1 ) return 1 ;

	std::cout << "dense_array construction (deep copy) from std::initializer_list and shape which can be value, std::initializer_list, Container, or Array" << '\n' ;
	std::cout << "glas3::dense_array<int> a6({1, 2, 3, 4, 5, 6}, {2, 3}) ;" << '\n' ;
	glas3::dense_array<double> a6({1, 2, 3, 4, 5, 6}, {2, 3}) ;
	if ( a6.size() != 6 || a6.ndof() != 6 || a6.shape().size() != 2 || a6.shape()[0] != 2 || a6.shape()[1] != 3 || a6.shape().shape().size() != 1 ) return 1 ;
	for ( i = 0; i < 2; ++i ) {
		for ( j = 0; j < 3; ++j ) {
			if ( a6[i + 2 * j] != i + 2 * j + 1 || a6({i, j}) != i + 2 * j + 1 ) return 1 ;
		}
	}

	std::cout << "dense_array construction (deep copy) from Contrainer std::vector and shape which can be value, std::initializer_list, Container, or Array" << '\n' ;
	std::cout << "glas3::dense_array<int> a9(std::vector<int>({1, 2, 3, 4, 5, 6}), std::vector<int>({2, 3})) ;" << '\n' ;
	glas3::dense_array<double> a9(std::vector<int>({1, 2, 3, 4, 5, 6}), std::vector<int>({2, 3})) ;
	if ( a9.size() != 6 || a9.ndof() != 6 || a9.shape().size() != 2 || a9.shape()[0] != 2 || a9.shape()[1] != 3 || a9.shape().shape().size() != 1 ) return 1 ;
	for ( i = 0; i < 2; ++i ) {
		for ( j = 0; j < 3; ++j ) {
			if ( a9[i + 2 * j] != i + 2 * j + 1 || a9({i, j}) != i + 2 * j + 1 ) return 1 ;
		}
	}

	std::cout << "dense_array construction (deep copy) from dense_array" << '\n' ;
	std::cout << "glas3::dense_array<double> a7( a6 ) ;" << '\n' ;
	glas3::dense_array<double> a7( a6 ) ;
	if ( a7.size() != a6.size() || a7.ndof() != a6.ndof() || a7.shape().size() != a6.shape().size() || a7.shape()[0] != a6.shape()[0] || a7.shape()[1] != a6.shape()[1] || a7.shape().shape().size() != 1 ) return 1 ;
    auto dum1 = a6[0] ;
	a6[0] = dum1 - 1 ;
	for ( i = 0; i < 2; ++i ) {
		for ( j = 0; j < 3; ++j ) {
			if ( i == 0 && j == 0 ) { if ( a7[i + 2 * j] != dum1 || a7({i, j}) != dum1 ) return 1 ; }
			else { if ( a7[i + 2 * j] != a6[i + 2 * j] || a7({i, j}) != a6({i, j}) ) return 1 ; }
		}
	}

	std::cout << '\n' ;

	return 0;
}
