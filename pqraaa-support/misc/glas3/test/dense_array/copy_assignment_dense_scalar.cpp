#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>

#include <iostream>
#include <vector>

int main() {

 	std::cout << '\n';

	std::cout << "copy_assignment of empty dense_scalar from dense_scalar, std::initializer_list, Container, or Array (deep copy)" << '\n' ;
	std::cout << "glas3::dense_scalar<double> s1 ;" << '\n' ;
	std::cout << "glas3::dense_scalar<double> s2 ;" << '\n' ;
	std::cout << "glas3::dense_vector<double> v1 ;" << '\n' ;
	std::cout << "s1 = s2 ;" << '\n' ;
	std::cout << "s1 = {} ;" << '\n' ;
	std::cout << "s1 = v1 ;" << '\n' ;
	std::cout << "s1 = std::vector<double>() ;" << '\n' ;
	glas3::dense_scalar<double> s1 ;
	glas3::dense_scalar<double> s2 ;
	glas3::dense_vector<double> v1 ;

	s1 = s2 ;
	if ( s1.size() != 0 || s1.ndof() != 0 || s1.shape().size() != 0 || s1.shape().shape().size() != 0 ) return 1 ;

	s1 = {} ;
	if ( s1.size() != 0 || s1.ndof() != 0 || s1.shape().size() != 0 || s1.shape().shape().size() != 0 ) return 1 ;

	s1 = v1 ;
	if ( s1.size() != 0 || s1.ndof() != 0 || s1.shape().size() != 0 || s1.shape().shape().size() != 0 ) return 1 ;

	s1 = std::vector<double>() ;
	if ( s1.size() != 0 || s1.ndof() != 0 || s1.shape().size() != 0 || s1.shape().shape().size() != 0 ) return 1 ;

	std::cout << "copy_assignment of dense_scalar from dense_scalar, value, std::initializer_list, Container, Array (deep copy)" << '\n' ;
	std::cout << "glas3::dense_scalar<double> s3 = 1 ;" << '\n' ;
	std::cout << "glas3::dense_scalar<double> s4 = 2 ;" << '\n' ;
	std::cout << "glas3::dense_vector<double> v2 = 3 ;" << '\n' ;
	std::cout << "s3 = s4 ;" << '\n' ;
	std::cout << "s3 = 5 ;" << '\n' ;
	std::cout << "s3 = {7} ;" << '\n' ;
	std::cout << "s3 = v2 ;" << '\n' ;
	std::cout << "s3 = std::vector<double>({5}) ;" << '\n' ;
	glas3::dense_scalar<double> s3 = 1 ;
	glas3::dense_scalar<double> s4 = 2 ;
	glas3::dense_vector<double> v2 = 3 ;

	s3 = s4 ;
	if ( s3.size() != 1 || s3.ndof() != 1 || s3.shape().size() != 0 || s3.shape().shape().size() != 0 ) return 1 ;
	auto dum = s3[0] ;
	s4[0] += 1 ;
	if ( s3[0] != dum ) return 1 ;

	s3 = 5 ;
	if ( s3.size() != 1 || s3.ndof() != 1 || s3.shape().size() != 0 || s3.shape().shape().size() != 0 ) return 1 ;
	if ( s3[0] != 5 ) return 1 ;

	s3 = {4} ;
	if ( s3.size() != 1 || s3.ndof() != 1 || s3.shape().size() != 0 || s3.shape().shape().size() != 0 ) return 1 ;
	if ( s3[0] != 4 ) return 1 ;

	s3 = v2 ;
	if ( s3.size() != 1 || s3.ndof() != 1 || s3.shape().size() != 0 || s3.shape().shape().size() != 0 ) return 1 ;
	dum = s3[0] ;
	v2[0] += 1 ;
	if ( s3[0] != dum ) return 1 ;

	s3 = std::vector<double>({5}) ;
	if ( s3.size() != 1 || s3.ndof() != 1 || s3.shape().size() != 0 || s3.shape().shape().size() != 0 ) return 1 ;
	if ( s3[0] != 5 ) return 1 ;

	std::cout << '\n' ;
	return 0;
}
