#include <glas3/array/dense_array/container/dense_vector.hpp>

#include <iostream>
#include <vector>

int main() {

 	std::cout << '\n';

	std::cout << "empty dense_vector construction" << '\n' ;
	std::cout << "glas3::dense_vector<double> v1 ;" << '\n' ;
	glas3::dense_vector<double> v1 ;
	if ( v1.size() != 0 || v1.ndof() != 0 || v1.shape().size() != 1 || v1.shape()[0] != 0 || v1.shape().shape().size() != 0 ) return 1 ;

	std::cout << "dense_vector construction (deep copy) from value" << '\n' ;
	std::cout << "glas3::dense_vector<int> v2(7) ;" << '\n' ;
	glas3::dense_vector<int> v2(7) ;
	if ( v2.size() != 1 || v2.ndof() != 1 || v2.shape().size() != 1 || v2.shape()[0] != 1 || v2.shape().shape().size() != 0 ) return 1 ;
	if ( v2[0] != 7 || v2({0}) != 7 ) return 1 ;

	std::cout << "dense_vector construction (deep copy) from std::initializer_list" << '\n' ;
	std::cout << "glas3::dense_vector<int> v3({3, 6, 1}) ;" << '\n' ;
	glas3::dense_vector<int> v3({3, 6, 1}) ;
	if ( v3.size() != 3 || v3.ndof() != 3 || v3.shape().size() != 1 || v3.shape()[0] != 3 || v3.shape().shape().size() != 0 ) return 1 ;
	if ( v3[0] != 3 || v3[1] != 6 || v3[2] != 1 || v3({0}) != 3 || v3({1}) != 6 || v3({2}) != 1 ) return 1 ;

	std::cout << "dense_vector construction (deep copy) from Container std::vector" << '\n' ;
	std::cout << "glas3::dense_vector<int> v7(std::vector<int>({3, 6, 1})) ;" << '\n' ;
	glas3::dense_vector<int> v7(std::vector<int>({3, 6, 1})) ;
	if ( v7.size() != 3 || v7.ndof() != 3 || v7.shape().size() != 1 || v7.shape()[0] != 3 || v7.shape().shape().size() != 0 ) return 1 ;
	if ( v7[0] != 3 || v7[1] != 6 || v7[2] != 1 || v7({0}) != 3 || v7({1}) != 6 || v7({2}) != 1 ) return 1 ;

	std::cout << "dense_vector construction (deep copy) from dense_vector" << '\n' ;
	std::cout << "glas3::dense_vector<int> v4(v3) ;" << '\n' ;
	glas3::dense_vector<int> v4(v3) ;
	if ( v4.size() != v3.size() || v4.ndof() !=  v3.ndof() || v4.shape().size() != v3.shape().size() || v4.shape()[0] != v3.shape()[0] || v4.shape().shape().size() != 0 ) return 1 ;
	auto dum1 = v3[0] ;
	v3[0] = dum1 - 1 ;
	if ( v4[0] != dum1 || v4[1] != v3[1] || v4[2] != v3[2] || v4({0}) != dum1 || v4({1}) != v3({1}) || v4({2}) != v3({2}) ) return 1 ;

	std::cout << "dense_vector construction (deep copy) from value and size" << '\n' ;
	std::cout << "glas3::dense_vector<int> v5(0, 2) ;" << '\n' ;
	glas3::dense_vector<int> v5(0, 2) ;
	if ( v5.size() != 2 || v5.ndof() != 2 || v5.shape().size() != 1 || v5.shape()[0] != 2 || v5.shape().shape().size() != 0 ) return 1 ;
	if ( v5[0] != 0 || v5[1] != 0 || v5({0}) != 0 || v5({1}) != 0 ) return 1 ;

	std::cout << "uninitialized dense_vector construction (deep copy) from size -> boost::make_shared_noinit apparently returns default initialized??" << '\n' ;
	std::cout << "glas3::dense_vector<int> v6(glas3::no_init(), 2) ;" << '\n' ;
	glas3::dense_vector<int> v6(glas3::no_init(), 2) ;
	if ( v6.size() != 2 || v6.ndof() != 2 || v6.shape().size() != 1 || v6.shape()[0] != 2 || v6.shape().shape().size() != 0 ) return 1 ;

	std::cout << '\n' ;
	return 0;
}
