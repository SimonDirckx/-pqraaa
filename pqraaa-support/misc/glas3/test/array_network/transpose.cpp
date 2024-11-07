#include <glas3/array/array_network/type/CP_tensor.hpp>
#include <glas3/array/array_network/type/Tucker_tensor.hpp>
#include <glas3/array/array_network/type/tensor_train.hpp>

#include <glas3/array/array_network/algorithm/randomize.hpp>
#include <glas3/array/array_network/algorithm/full.hpp>
#include <glas3/array/array_network/algorithm/transpose.hpp>

#include <glas3/dense_array.hpp>

#include <random>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <iostream>

int main() {

 	std::cout << '\n';

 	typedef typename boost::uniform_real<double>                        distribution_type ;
 	typedef typename boost::lagged_fibonacci607                         engine_type ;
 	typedef boost::variate_generator< engine_type, distribution_type >  generator_type ;
 	std::random_device rd ;
 	generator_type g( engine_type( rd() ), distribution_type( -1, 1 ) ) ;

    std::cout << "CP tensor construction (deep copy) from value, shape and rank" << '\n' ;
    std::cout << "glas3::CP_tensor<glas3::dense_array<double>> t1( 0, {3, 4, 5}, 3 ) ;" << '\n' << '\n' ;
    glas3::dense_vector<std::ptrdiff_t> s1( {3, 4, 5} ) ;
    std::ptrdiff_t r1( 3 ) ;
    glas3::CP_tensor<glas3::dense_array<double>> t1( 0, s1, r1 ) ;
    glas3::randomize( t1, g ) ;
    std::cout << "transpose operation (shallow copy)" << '\n' ;
    std::cout << "auto t1_p = glas3::transpose( t1) ;" << '\n' << '\n' ;
    auto t1_p = glas3::transpose( t1 ) ;
    glas3::randomize( t1, g ) ;
    auto t1_full_p = glas3::transpose( glas3::full( t1 ) );

    if ( glas3::norm_inf( glas3::full( t1_p ) - t1_full_p ) > 1e-14 ) return 1 ;

    std::cout << "tensor_train construction (deep copy) from value, shape and rank" << '\n' ;
    std::cout << "glas3::tensor_train<glas3::dense_array<double>> t2(0, {5, 4, 3}, {4, 3}) ;" << '\n' << '\n' ;
    glas3::dense_vector<std::ptrdiff_t> s2( {5, 4, 3} ) ;
    glas3::dense_vector<std::ptrdiff_t> r2( {4, 3} ) ;
    glas3::tensor_train<glas3::dense_array<double>> t2(0, s2, r2) ;
    glas3::randomize( t2, g ) ;
    std::cout << "transpose operation (shallow copy)" << '\n' ;
    std::cout << "auto t2_p = glas3::transpose( t2 ) ;" << '\n' ;
    auto t2_p = glas3::transpose( t2 ) ;
    glas3::randomize( t2, g ) ;
    auto t2_full_p = glas3::transpose( glas3::full( t2 ) );

    if ( glas3::norm_inf( glas3::full( t2_p ) - t2_full_p ) > 1e-14 ) return 1 ;

    std::cout << '\n' ;

	return 0;

}
