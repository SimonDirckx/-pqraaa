#include <glas3/array/dense_array/container/dense_matrix.hpp>

#include <iostream>
#include <vector>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j ;

	std::cout << "empty dense_matrix construction" << '\n' ;
	std::cout << "glas3::dense_matrix<double> m1 ;" << '\n' ;
	glas3::dense_matrix<double> m1 ;
	if ( m1.size() != 0 || m1.ndof() != 0 || m1.shape().size() != 2 || m1.shape()[0] != 0 || m1.shape()[1] != 0 || m1.shape().shape().size() != 1 ) return 1 ;

	std::cout << "dense_matrix construction (deep copy) from value" << '\n' ;
	std::cout << "glas3::dense_matrix<int> m2(7) ;" << '\n' ;
	glas3::dense_matrix<int> m2(7) ;
	if ( m2.size() != 1 || m2.ndof() != 1 || m2.shape().size() != 2 || m2.shape()[0] != 1 || m2.shape()[1] != 1 || m2.shape().shape().size() != 1 ) return 1 ;
	if ( m2[0] != 7 || m2({0, 0}) != 7 ) return 1 ;

	std::cout << "dense_matrix construction (deep copy) from std::initializer_list" << '\n' ;
	std::cout << "glas3::dense_matrix<int> m3({1, 2}) ;" << '\n' ;
	glas3::dense_matrix<int> m3({1, 2}) ;
	if ( m3.size() != 2 || m3.ndof() != 2 || m3.shape().size() != 2 || m3.shape()[0] != 2 || m3.shape()[1] != 1 || m3.shape().shape().size() != 1 ) return 1 ;
    if ( m3[0] != 1 || m3[1] != 2 || m3({0, 0}) != 1 || m3({1, 0}) != 2 ) return 1 ;

	std::cout << "dense_matrix construction (deep copy) from Container std::vector" << '\n' ;
	std::cout << "glas3::dense_matrix<int> m8(std::vector<double>({1, 2})) ;" << '\n' ;
	glas3::dense_matrix<int> m8(std::vector<double>({1, 2})) ;
	if ( m8.size() != 2 || m8.ndof() != 2 || m8.shape().size() != 2 || m8.shape()[0] != 2 || m8.shape()[1] != 1 || m8.shape().shape().size() != 1 ) return 1 ;
    if ( m8[0] != 1 || m8[1] != 2 || m8({0, 0}) != 1 || m8({1, 0}) != 2 ) return 1 ;

    std::cout << "dense_matrix construction (deep copy) from value and shape which can be value, std::initializer_list, Container, or Array" << '\n' ;
    std::cout << "glas3::dense_matrix<int> m4(0, {2, 3}) ;" << '\n' ;
    glas3::dense_matrix<int> m4(0, {2, 3}) ;
    if ( m4.size() != 6 || m4.ndof() != 6 || m4.shape().size() != 2 || m4.shape()[0] != 2 || m4.shape()[1] != 3 || m4.shape().shape().size() != 1 ) return 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		if ( m4[i + 2 * j] != 0 || m4({i, j}) != 0 ) return 1 ;
    	}
    }

    std::cout << "uninitialized dense_matrix construction (deep copy) from shape -> boost::make_shared_noinit apparently returns default initialized??" << '\n' ;
    std::cout << "glas3::dense_matrix<int> m5(glas3::no_init(), {2, 3}) ;" << '\n' ;
    glas3::dense_matrix<int> m5(glas3::no_init(), {2, 3}) ;
    if ( m5.size() != 6 || m5.ndof() != 6 || m5.shape().size() != 2 || m5.shape()[0] != 2 || m5.shape()[1] != 3 || m5.shape().shape().size() != 1 ) return 1 ;

	std::cout << "dense_matrix construction (deep copy) from std::initializer_list and shape which can be value, std::initializer_list, Container, or Array" << '\n' ;
	std::cout << "glas3::dense_matrix<int> m6({1, 2, 3, 4, 5, 6}, {2, 3}) ;" << '\n' ;
	glas3::dense_matrix<double> m6({1, 2, 3, 4, 5, 6}, {2, 3}) ;
	if ( m6.size() != 6 || m6.ndof() != 6 || m6.shape().size() != 2 || m6.shape()[0] != 2 || m6.shape()[1] != 3 || m6.shape().shape().size() != 1 ) return 1 ;
	for ( i = 0; i < 2; ++i ) {
		for ( j = 0; j < 3; ++j ) {
			if ( m6[i + 2 * j] != i + 2 * j + 1 || m6({i, j}) != i + 2 * j + 1 ) return 1 ;
		}
	}

	std::cout << "dense_matrix construction (deep copy) from Contrainer std::vector and shape which can be value, std::initializer_list, Container, or Array" << '\n' ;
	std::cout << "glas3::dense_matrix<int> m9(std::vector<double>({1, 2, 3, 4, 5, 6}), std::vector<double>({2, 3})) ;" << '\n' ;
	glas3::dense_matrix<double> m9(std::vector<double>({1, 2, 3, 4, 5, 6}), std::vector<double>({2, 3})) ;
	if ( m9.size() != 6 || m9.ndof() != 6 || m9.shape().size() != 2 || m9.shape()[0] != 2 || m9.shape()[1] != 3 || m9.shape().shape().size() != 1 ) return 1 ;
	for ( i = 0; i < 2; ++i ) {
		for ( j = 0; j < 3; ++j ) {
			if ( m9[i + 2 * j] != i + 2 * j + 1 || m9({i, j}) != i + 2 * j + 1 ) return 1 ;
		}
	}

	std::cout << "dense_matrix construction (deep copy) from dense_matrix" << '\n' ;
	std::cout << "glas3::dense_matrix<double> m7(m6) ;" << '\n' ;
	glas3::dense_matrix<double> m7(m6) ;
	if ( m7.size() != m6.size() || m7.ndof() != m6.ndof() || m7.shape().size() != m6.shape().size() || m7.shape()[0] != m6.shape()[0] || m7.shape()[1] != m6.shape()[1] || m7.shape().shape().size() != 1 ) return 1 ;
    auto dum1 = m6[0] ;
	m6[0] = dum1 - 1 ;
	for ( i = 0; i < 2; ++i ) {
		for ( j = 0; j < 3; ++j ) {
			if ( i == 0 && j == 0 ) { if ( m7[i + 2 * j] != dum1 || m7({i, j}) != dum1 ) return 1 ; }
			else { if ( m7[i + 2 * j] != m6[i + 2 * j] || m7({i, j}) != m6({i, j}) ) return 1 ; }
		}
	}

	std::cout << '\n' ;

	return 0;
}
