#include <glas3/glas3.hpp>

#include <iostream>
#include <cmath>

int main() {

 	std::cout << '\n';

 	glas3::dense_array<double> a1({-3, 0, -2, 0, -1, 9}, {2, 3, 1}) ;
 	glas3::dense_array<std::complex<double>> a2({-3, 0, std::complex<double>(0, -2), 0, std::complex<double>(0, -1), 9}, {3, 2}) ;

 	std::cout << "glas3::dense_array<double> a1({-3, 0, -2, 0, -1, 9}, {2, 3, 1}) ;" << '\n' << a1 << '\n' ;
 	std::cout << "glas3::dense_array<std::complex<double>> a2({-3, 0, std::complex<double>(0, -2), 0, std::complex<double>(0, -1), 9}, {3, 2}) ;" << '\n' << a2 << '\n' ;
 	std::cout << "glas3::norm(a1, 0) ;" << '\n' << glas3::norm(a1, 0) << '\n' << '\n' ;
 	std::cout << "glas3::norm_2(a2) ;" << '\n' << glas3::norm_2(a2) << '\n' << '\n' ;
 	std::cout << "glas3::norm(a2, INFINITY) ;" << '\n' << glas3::norm(a2, INFINITY) << '\n' << '\n' ;
 	std::cout << "glas3::inner_prod(a1, a2) ;" << '\n' << glas3::inner_prod(a1, a2) << '\n' << '\n' ;
 	std::cout << "other operations:" << '\n' << "all, any, assign, fill, find, ndims, ndof, prod, shape, size, sort, sum" << '\n' << '\n' ;

	return 0;
}
