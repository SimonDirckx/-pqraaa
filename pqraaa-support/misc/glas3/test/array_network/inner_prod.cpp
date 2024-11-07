#include <glas3/array/array_network/type/CP_tensor.hpp>
#include <glas3/array/array_network/type/Tucker_tensor.hpp>
#include <glas3/array/array_network/type/tensor_train.hpp>

#include <glas3/array/array_network/algorithm/randomize.hpp>
#include <glas3/array/array_network/algorithm/full.hpp>
#include <glas3/array/array_network/algorithm/inner_prod.hpp>

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
    auto t1_full = glas3::full( t1 ) ;

    std::cout << "tensor_train construction (deep copy) from value, shape and rank" << '\n' ;
    std::cout << "glas3::tensor_train<glas3::dense_array<double>> t2(0, {5, 4, 3}, {4, 3}) ;" << '\n' ;
    glas3::dense_vector<std::ptrdiff_t> s2( {3, 4, 5} ) ;
    glas3::dense_vector<std::ptrdiff_t> r2( {4, 3} ) ;
    glas3::tensor_train<glas3::dense_array<double>> t2(0, s2, r2) ;
    glas3::randomize( t2, g ) ;
    auto t2_full = glas3::full( t2 ) ;

    std::cout << "Tucker tensor construction (deep copy) from value, shape and rank -> boost::make_shared_noinit apparently returns default initialized??" << '\n' ;
    std::cout << "glas3::Tucker_tensor<glas3::dense_array<double>> t3(1, {4, 3, 5}, {4, 3, 2}) ;" << '\n' << '\n' ;
    glas3::dense_vector<std::ptrdiff_t> s3( {3, 4, 5} ) ;
    glas3::dense_vector<std::ptrdiff_t> r3( {4, 3, 2} ) ;
    glas3::Tucker_tensor<glas3::dense_array<double>> t3(0, s3, r3) ;
    glas3::randomize( t3, g ) ;
    auto t3_full = glas3::full( t3 ) ;

	std::cout << "inner product operation" << '\n' ;
	std::cout << "auto t2_t3_i = glas3::inner_prod( t2, t3 ) ;" << '\n' ;
	std::cout << "auto t2_t1_i = glas3::inner_prod( t2, t1 ) ;" << '\n' ;
	std::cout << "auto t3_t1_i = glas3::inner_prod( t3, t1 ) ;" << '\n' ;
	auto t2_t3_i = glas3::inner_prod( t2, t3 ) ;
	auto t2_t1_i = glas3::inner_prod( t2, t1 ) ;
	auto t3_t1_i = glas3::inner_prod( t3, t1 ) ;
	auto t2_t3_i_ref = glas3::inner_prod( t2_full, t3_full ) ;
	auto t2_t1_i_ref = glas3::inner_prod( t2_full, t1_full ) ;
	auto t3_t1_i_ref = glas3::inner_prod( t3_full, t1_full ) ;

	std::cout << "std::abs( t2_t3_i - t2_t3_i_ref ): " << std::abs( t2_t3_i - t2_t3_i_ref ) << '\n' ;
	std::cout << "std::abs( t2_t1_i - t2_t1_i_ref ): " << std::abs( t2_t1_i - t2_t1_i_ref ) << '\n' ;
	std::cout << "std::abs( t3_t1_i - t3_t1_i_ref ): " << std::abs( t3_t1_i - t3_t1_i_ref ) << '\n' ;
    if ( std::abs( t2_t3_i - t2_t3_i_ref ) > 1e-14 ) return 1 ;
    if ( std::abs( t2_t1_i - t2_t1_i_ref ) > 1e-14 ) return 1 ;
    if ( std::abs( t3_t1_i - t3_t1_i_ref ) > 1e-14 ) return 1 ;

	std::cout << '\n' ;

	return 0;

}
