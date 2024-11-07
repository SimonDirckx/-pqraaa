#include <glas3/glas3.hpp>

#include <iostream>
#include <cmath>

int main() {

 	std::cout << '\n';

 	glas3::dense_array<double> a(glas3::no_init(), {2, 3, 2}) ;
 	a = glas3::range<>(0, 12) ;

 	std::cout << "glas3::dense_array<double> a(glas3::no_init(), {2, 3, 2}) ;" << '\n' << "a = glas3::range<>(0, 12) ;" << '\n' << a << '\n' ;
 	std::cout << "a.size() ;" << '\n' << a.size() << '\n' << '\n' ;
 	std::cout << "a.shape() ;" << '\n' << a.shape() << '\n' ;
 	std::cout << "a.ndof() ;" << '\n' << a.ndof() << '\n' << '\n' ;
 	std::cout << "a.shape().size() ;" << '\n' << a.shape().size() << '\n' << '\n' ;
 	std::cout << "a[2] ;" << '\n' << a[2] << '\n' << '\n' ;
 	std::cout << "a({1, 2, 1}) ;" << '\n' << a({1, 2, 1}) << '\n' << '\n' ;
 	auto b = a.shallow_copy() ;
 	a[0] = -1 ;
 	std::cout << "auto b = a.shallow_copy() ;" << '\n' ;
 	std::cout << "a[0] = -1 ;" << '\n' ;
 	std::cout << "b ;" << '\n' << b << '\n' ;

	return 0;
}
