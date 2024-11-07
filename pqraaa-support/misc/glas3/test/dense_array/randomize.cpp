#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/iostream.hpp>
#include <glas3/array/dense_array/algorithm/randomize.hpp>

#include <random>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/lagged_fibonacci.hpp>

#include <complex>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::cout << "randomize operation" << '\n' << '\n' ;

	std::cout << "std::random_device rd ;" << '\n' ;
 	std::cout << "boost::variate_generator<..., ...> g( boost::lagged_fibonacci607( rd() ), boost::uniform_real<double>( -1, 1 ) ) ;" << '\n' << '\n' ;
	typedef typename boost::uniform_real<double>                        distribution_type ;
	typedef typename boost::lagged_fibonacci607                         engine_type ;
	typedef boost::variate_generator< engine_type, distribution_type >  generator_type ;
 	std::random_device rd ;
	generator_type g( engine_type( rd() ), distribution_type( -1, 1 ) ) ;

	std::cout << "glas3::dense_array<double> a(glas3::no_init(), {2, 3}) ;" << '\n' ;
 	std::cout << "glas3::randomize(a, g) ;" << '\n' ;
 	glas3::dense_array<double> a(glas3::no_init(), {2, 3}) ;
 	glas3::randomize(a, g) ;
    std::cout << a << '\n' ;

 	std::cout << "glas3::dense_array<std::complex<double>> b(glas3::no_init(), {2, 3}) ;" << '\n' ;
 	std::cout << "glas3::randomize(b, g) ;" << '\n' ;
    glas3::dense_array<std::complex<double>> b(glas3::no_init(), {2, 3}) ;
    glas3::randomize(b, g) ;
    std::cout << b ;

    std::cout << '\n' ;

	return 0;

}
