#include <glas3/array/dense_array/type/range.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

	std::cout << "range construction" << '\n' ;
	std::cout << "glas3::range<> r1(10, 10) ;" << '\n' ;
	glas3::range<> r1(10, 10) ;
	if ( r1.size() != 0 || r1.ndof() != 3 || r1.shape().size() != 1 || r1.shape()[0] != 0 || r1.shape().shape().size() != 0 ) return 1 ;

	std::cout << "range construction" << '\n' ;
	std::cout << "glas3::range<> r2(10, 20) ;" << '\n' ;
	glas3::range<> r2(10, 13) ;
	if ( r2.size() != 3 || r2.ndof() != 3 || r2.shape().size() != 1 || r2.shape()[0] != 3 || r2.shape().shape().size() != 0 ) return 1 ;
	if ( r2[0] != 10 || r2[1] != 11 || r2[2] != 12 || r2({0}) != 10 || r2({1}) != 11 || r2({2}) != 12 ) return 1 ;

	std::cout << "empty range construction" << '\n' ;
	std::cout << "glas3::range<> r3(10, 1) ;" << '\n' ;
	glas3::range<> r3(10, 1) ;
	if ( r3.size() != 0 || r3.ndof() != 3 || r3.shape().size() != 1 || r3.shape()[0] != 0 || r3.shape().shape().size() != 0 ) return 1 ;

	std::cout << "range construction" << '\n' ;
	std::cout << "glas3::range<> r4(0, 4, 2) ;" << '\n' ;
	glas3::range<> r4(0, 4, 2) ;
	if ( r4.size() != 2 || r4.ndof() != 3 || r4.shape().size() != 1 || r4.shape()[0] != 2 || r4.shape().shape().size() != 0 ) return 1 ;
	if ( r4[0] != 0 || r4[1] != 2 || r4({0}) != 0 || r4({1}) != 2 ) return 1 ;

	std::cout << "range construction" << '\n' ;
	std::cout << "glas3::range<> r5(0, 5, 2) ;" << '\n' ;
	glas3::range<> r5(0, 5, 2) ;
	if ( r5.size() != 3 || r5.ndof() != 3 || r5.shape().size() != 1 || r5.shape()[0] != 3 || r5.shape().shape().size() != 0 ) return 1 ;
	if ( r5[0] != 0 || r5[1] != 2 || r5[2] != 4 || r5({0}) != 0 || r5({1}) != 2 || r5({2}) != 4 ) return 1 ;

	std::cout << "range construction" << '\n' ;
	std::cout << "glas3::range<> r6(0, -4, -2) ;" << '\n' ;
	glas3::range<> r6(0, -4, -2) ;
	if ( r6.size() != 2 || r6.ndof() != 3 || r6.shape().size() != 1 || r6.shape()[0] != 2 || r6.shape().shape().size() != 0 ) return 1 ;
	if ( r6[0] != 0 || r6[1] != -2 || r6({0}) != 0 || r6({1}) != -2 ) return 1 ;

	std::cout << "range construction" << '\n' ;
	std::cout << "glas3::range<> r7(0, -5, -2) ;" << '\n' ;
	glas3::range<> r7(0, -5, -2) ;
	if ( r7.size() != 3 || r7.ndof() != 3 || r7.shape().size() != 1 || r7.shape()[0] != 3 || r7.shape().shape().size() != 0 ) return 1 ;
	if ( r7[0] != 0 || r7[1] != -2 || r7[2] != -4 || r7({0}) != 0 || r7({1}) != -2 || r7({2}) != -4 ) return 1 ;

	std::cout << "empty range construction" << '\n' ;
	std::cout << "glas3::range<> r8(0, -5, 2) ;" << '\n' ;
	glas3::range<> r8(0, -5, 2) ;
	if ( r8.size() != 0 || r8.ndof() != 3 || r8.shape().size() != 1 || r8.shape()[0] != 0 || r8.shape().shape().size() != 0 ) return 1 ;

	std::cout << '\n' ;
	return 0;
}
