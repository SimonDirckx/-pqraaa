#include <glas3/array/type/empty_array.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::cout << "empty_array construction" << '\n' ;
	std::cout << "glas3::empty_array e1 ;" << '\n' ;
	glas3::empty_array e1 ;
	if ( e1.size() != 0 || e1.ndof() != 0 || e1.shape().size() != 0 || e1.shape().shape().size() != 0 ) return 1 ;

	std::cout << "empty_array construction (deep copy) from std::initializer_list" << '\n' ;
	std::cout << "glas3::empty_array e2 ;" << '\n' ;
	glas3::empty_array e2 = {} ;
	if ( e2.size() != 0 || e2.ndof() != 0 || e2.shape().size() != 0 || e2.shape().shape().size() != 0 ) return 1 ;

	std::cout << '\n' ;
	return 0;

}
