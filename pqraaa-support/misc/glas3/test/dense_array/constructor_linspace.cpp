#include <glas3/array/dense_array/type/linspace.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

	std::cout << "linspace construction" << '\n' ;
	std::cout << "glas3::linspace<double> l1(0, 4, 3) ;" << '\n' ;
	glas3::linspace<double> l1(0, 4, 3) ;
	if ( l1.size() != 3 || l1.ndof() != 3 || l1.shape().size() != 1 || l1.shape()[0] != 3 || l1.shape().shape().size() != 0 ) return 1 ;
	if ( l1[0] != 0 || l1[1] != 2 || l1[2] != 4 || l1({0}) != 0 || l1({1}) != 2 || l1({2}) != 4 ) return 1 ;

	std::cout << "empty linspace construction" << '\n' ;
	std::cout << "glas3::linspace<double> l2(0, 4, 0) ;" << '\n' ;
	glas3::linspace<double> l2(0, 4, 0) ;
	if ( l2.size() != 0 || l2.ndof() != 3 || l2.shape().size() != 1 || l2.shape()[0] != 0 || l2.shape().shape().size() != 0 ) return 1 ;

	std::cout << "linspace construction" << '\n' ;
	std::cout << "glas3::linspace<double> l3(0, -4, 3) ;" << '\n' ;
	glas3::linspace<double> l3(0, -4, 3) ;
	if ( l3.size() != 3 || l3.ndof() != 3 || l3.shape().size() != 1 || l3.shape()[0] != 3 || l3.shape().shape().size() != 0 ) return 1 ;
	if ( l3[0] != 0 || l3[1] != -2 || l3[2] != -4 || l3({0}) != 0 || l3({1}) != -2 || l3({2}) != -4 ) return 1 ;

	std::cout << '\n' ;
	return 0;
}
