#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/type/scalar_wrapper.hpp>

#include <iostream>
#include <vector>

int main() {

 	std::cout << '\n';

 	std::cout << "empty scalar_wrapper construction" << '\n' ;
	std::cout << "glas3::scalar_wrapper<double> w1 ;" << '\n' ;
	glas3::scalar_wrapper<double> w1 ;
	if ( w1.size() != 0 || w1.ndof() != 0 || w1.shape().size() != 0 || w1.shape().shape().size() != 0 ) return 1 ;

	std::cout << "scalar_wrapper construction (deep copy) from value" << '\n' ;
	std::cout << "glas3::scalar_wrapper<double> w2( 4 ) ;" << '\n' ;
	glas3::scalar_wrapper<double> w2( 4 ) ;
	if ( w2.size() != 1 || w2.ndof() != 1 || w2.shape().size() != 0 || w2.shape().shape().size() != 0 ) return 1 ;
	if ( w2[0] != 4 ) return 1 ;

	std::cout << "scalar_wrapper construction (deep copy) from std::initializer_list" << '\n' ;
	std::cout << "glas3::scalar_wrapper<double> w3( {4} ) ;" << '\n' ;
	glas3::scalar_wrapper<double> w3( {4} ) ;
	if ( w3.size() != 1 || w3.ndof() != 1 || w3.shape().size() != 0 || w3.shape().shape().size() != 0 ) return 1 ;
	if ( w3[0] != 4 ) return 1 ;

	std::cout << "scalar_wrapper construction (shallow copy) from Array" << '\n' ;
	std::cout << "glas3::scalar_wrapper<double> w6( glas3::dense_scalar<double>( 5 ) ) ;" << '\n' ;
	glas3::dense_scalar<double> s( 5 ) ;
	glas3::scalar_wrapper<double> w6( s ) ;
	if ( w6.size() != 1 || w6.ndof() != 1 || w6.shape().size() != 0 || w6.shape().shape().size() != 0 ) return 1 ;
	s[0] += 1 ;
	if ( w6[0] != s[0] ) return 1 ;

//	std::cout << "scalar_wrapper construction (deep copy) from Container std::vector" << '\n' ;
//	std::cout << "glas3::scalar_wrapper<double> w4( std::vector<double>( 1, 4 ) ) ;" << '\n' ;
//	glas3::scalar_wrapper<double> w4( std::vector<double>( 1, 4 ) ) ;
//	if ( w4.size() != 1 || w4.ndof() != 1 || w4.shape().size() != 0 || w4.shape().shape().size() != 0 ) return 1 ;
//	if ( w4[0] != 4 ) return 1 ;
//
//	std::cout << "scalar_wrapper construction (move) from RandomAccessContainer std::vector" << '\n' ;
//	std::cout << "glas3::scalar_wrapper<double> w5( std::move( std::vector<double>( 1, 4 ) ) ) ;" << '\n' ;
//	glas3::scalar_wrapper<double> w5( std::move( std::vector<double>( 1, 4 ) ) ) ;
//	if ( w5.size() != 1 || w5.ndof() != 1 || w5.shape().size() != 0 || w5.shape().shape().size() != 0 ) return 1 ;
//	if ( w5[0] != 4 ) return 1 ;

	std::cout << '\n' ;
	return 0;

}
