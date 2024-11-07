#include <glas3/array/type/empty_array.hpp>
#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_matrix.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <iostream>

int main() {

	std::ptrdiff_t i, j ;

 	std::cout << '\n';

	std::cout << "empty_array construction (deep copy) from Array" << '\n' ;
	std::cout << "glas3::empty_array e3 = glas3::dense_array<double>() ;" << '\n' ;
	auto dum2 = glas3::dense_array<double>() ;
	glas3::empty_array e3 = dum2 ;
	if ( e3.size() != dum2.size() || e3.ndof() != dum2.ndof() || e3.shape().size() != dum2.shape().size() || e3.shape().shape().size() != 0 ) return 1 ;

	std::cout << "dense_scalar construction (deep copy) from Array" << '\n' ;
	std::cout << "glas3::dense_scalar<int> s6 = glas3::dense_array<double>(-5.3) ;" << '\n' ;
	auto dum3 = glas3::dense_array<double>(-5.3) ;
	glas3::dense_scalar<double> s6 = dum3 ;
	if ( s6.size() != dum3.size() || s6.ndof() != dum3.ndof() || s6.shape().size() != dum3.shape().size() || s6.shape().shape().size() != 0 ) return 1 ;
	auto dum1 = dum3[0] ;
	dum3[0] = dum1 - 1 ;
	if ( s6[0] != dum1 ) return 1 ;

	std::cout << "dense_vector construction (deep copy) from Array" << '\n' ;
	std::cout << "glas3::dense_vector<int> v7 = glas3::dense_array<double>({0, 1, 2, 5}, {2, 2}) ;" << '\n' ;
	auto dum4 = glas3::dense_array<double>({0, 1, 2, 5}, {2, 2}) ;
	glas3::dense_vector<int> v7 = dum4 ;
	if ( v7.size() != dum4.size() || v7.ndof() != dum4.ndof() || v7.shape().size() != 1 || v7.shape()[0] != dum4.size() || v7.shape().shape().size() != 0 ) return 1 ;
	dum1 = dum4[0] ;
	dum4[0] = dum1 - 1 ;
	for ( i = 0; i < 4; ++i ) {
		if ( i == 0 ) { if ( v7[i] != dum1 || v7({i}) != dum1 ) return 1 ; }
		else { if ( v7[i] != dum4[i] || v7({i}) != dum4[i] ) return 1 ; }
	}

	std::cout << "dense_matrix construction (deep copy) from Array" << '\n' ;
	std::cout << "glas3::dense_matrix<double> m8 = glas3::dense_array<int>({0, 1, 2, 5}, {2, 2}) ;" << '\n' ;
	auto dum6 = glas3::dense_array<int>({0, 1, 2, 5}, {2, 2}) ;
	glas3::dense_matrix<double> m8 = dum6 ;
	if ( m8.size() != dum6.size() || m8.ndof() != dum6.ndof() || m8.shape().size() != dum6.shape().size() || m8.shape()[0] != dum6.shape()[0] || m8.shape()[1] != dum6.shape()[1] || m8.shape().shape().size() != 1 ) return 1 ;
	dum1 = dum6[0] ;
	dum6[0] = dum1 - 1 ;
	for ( i = 0; i < 2; ++i ) {
		for ( j = 0; j < 2; ++j ) {
			if ( i == 0 && j == 0 ) { if ( m8[i + 2 * j] != dum1 || m8({i, j}) != dum1 ) return 1 ; }
			else { if ( m8[i + 2 * j] != dum6[i + 2 * j] || m8({i, j}) != dum6({i, j}) ) return 1 ; }
		}
	}

	std::cout << "dense_array construction (deep copy) from Array" << '\n' ;
	std::cout << "glas3::dense_array<double> a8 = glas3::dense_scalar<int>(5) ;" << '\n' ;
	auto dum5 = glas3::dense_scalar<int>(5) ;
	glas3::dense_array<double> a8 = dum5 ;
	if ( a8.size() != dum5.size() || a8.ndof() != dum5.ndof() || a8.shape().size() != dum5.shape().size() || a8.shape().shape().size() != 1 ) return 1 ;
    dum1 = dum5[0] ;
    dum5[0] = dum1 - 1 ;
    if ( a8[0] != dum1 || a8({}) != dum1 ) return 1 ;

	std::cout << '\n' ;
	return 0;

}
