#include <glas3/array/dense_array/container/dense_scalar.hpp>

#include <iostream>
#include <vector>

int main() {

 	std::cout << '\n';

	std::cout << "empty dense_scalar construction" << '\n' ;
	std::cout << "glas3::dense_scalar<double> s1 ;" << '\n' ;
	glas3::dense_scalar<double> s1 ;
	if ( s1.size() != 0 || s1.ndof() != 0 || s1.shape().size() != 0 || s1.shape().shape().size() != 0 ) return 1 ;

    std::cout << "empty dense_scalar construction (deep copy) from std::initializer_list" << '\n' ;
    std::cout << "glas3::dense_scalar<double> s2({}) ;" << '\n' ;
	glas3::dense_scalar<double> s2 = {} ;
	if ( s2.size() != 0 || s2.ndof() != 0 || s2.shape().size() != 0 || s2.shape().shape().size() != 0 ) return 1 ;

    std::cout << "dense_scalar construction (deep copy) from value" << '\n' ;
    std::cout << "glas3::dense_scalar<double> s3(3.2) ;" << '\n' ;
	glas3::dense_scalar<double> s3(3.2) ;
	if ( s3.size() != 1 || s3.ndof() != 1 || s3.shape().size() != 0 || s3.shape().shape().size() != 0 ) return 1 ;
	if ( s3[0] != 3.2 ) return 1 ;

	std::cout << "dense_scalar construction (deep copy) from std::initializer_list" << '\n' ;
	std::cout << "glas3::dense_scalar<double> s4({-7.8}) ;" << '\n' ;
	glas3::dense_scalar<double> s4({-7.8}) ;
	if ( s4.size() != 1 || s4.ndof() != 1 || s4.shape().size() != 0 || s4.shape().shape().size() != 0 ) return 1 ;
	if ( s4[0] != -7.8 ) return 1 ;

	std::cout << "dense_scalar construction (deep copy) from Container std::vector" << '\n' ;
	std::cout << "glas3::dense_scalar<double> s6(std::vector<double>(1, -7.8)) ;" << '\n' ;
	glas3::dense_scalar<double> s6(std::vector<double>(1, -7.8)) ;
	if ( s6.size() != 1 || s6.ndof() != 1 || s6.shape().size() != 0 || s6.shape().shape().size() != 0 ) return 1 ;
	if ( s6[0] != -7.8 ) return 1 ;

	std::cout << "dense_scalar construction (deep copy) from dense_scalar" << '\n' ;
	std::cout << "glas3::dense_scalar<int> s5(s3) ;" << '\n' ;
	glas3::dense_scalar<double> s5(s3) ;
	if ( s5.size() != 1 || s5.ndof() != 1 || s5.shape().size() != 0 || s5.shape().shape().size() != 0 ) return 1 ;
	auto dum1 = s3[0] ;
	s3[0] = dum1 - 1 ;
	if ( s5[0] != dum1 ) return 1 ;

	std::cout << '\n' ;
	return 0;
}
