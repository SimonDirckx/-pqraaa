#include <glas3/glas3.hpp>

#include <iostream>
#include <cmath>

int main() {

 	std::cout << '\n';

 	glas3::empty_array e ;

 	glas3::dense_scalar<double> s = 3 ;
 	glas3::dense_vector<std::complex<double>> v = {-3, 0, std::complex<double>(3, -2)} ;
 	glas3::dense_array<double> a1(0, {2, 3}) ;
 	glas3::dense_array<double> a2(glas3::no_init(), {2, 3}) ;
 	glas3::dense_array<double> a3({1, 2, 3, 4, 5, 6}, {2, 3}) ;
 	glas3::dense_array<std::complex<double>> a4 = v ;

 	glas3::constant_array<double> c(0, {2, 3}) ;
 	glas3::eye<double> i({3, 3}) ;
 	glas3::linspace<double> l(0, -4, 3) ;
 	glas3::range<> r1(0, 4) ;
 	glas3::range<> r2(4, -1, -1) ;

 	std::cout << "glas3::empty_array e ;" << '\n' << e << '\n' ;
 	std::cout << "glas3::dense_scalar<double> s = 3 ;" << '\n' << s << '\n' ;
 	std::cout << "glas3::dense_vector<std::complex<double>> v = {-3, 0, std::complex<double>(3, -2)} ;" << '\n' << v << '\n' ;
 	std::cout << "glas3::dense_array<double> a1(0, {2, 3}) ;" << '\n' << a1 << '\n' ;
 	std::cout << "glas3::dense_array<double> a2(glas3::no_init(), {2, 3}) ;" << '\n' << a2 << '\n' ;
 	std::cout << "glas3::dense_array<double> a3({1, 2, 3, 4, 5, 6}, {2, 3}) ;" << '\n' << a3 << '\n' ;
 	std::cout << "glas3::dense_array<std::complex<double>> a4 = v ;" << '\n' << a4 << '\n' ;
 	std::cout << "glas3::constant_array<double> c(0, {2, 3}) ;" << '\n' << c << '\n' ;
 	std::cout << "glas3::eye<double> i({3, 3}) ;" << '\n' << i << '\n' ;
 	std::cout << "glas3::linspace<double> l(0, -4, 3) ;" << '\n' << l << '\n' ;
 	std::cout << "glas3::range<> r1(0, 4) ;" << '\n' << r1 << '\n' ;
 	std::cout << "glas3::range<> r2(4, -1, -1) ;" << '\n' << r2 << '\n' ;

	return 0;
}
