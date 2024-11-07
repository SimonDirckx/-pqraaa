#include <glas3/array/array_network/type/CP_tensor.hpp>
#include <glas3/array/array_network/type/Tucker_tensor.hpp>
#include <glas3/array/array_network/type/tensor_train.hpp>

#include <glas3/array/array_network/algorithm/randomize.hpp>
#include <glas3/array/array_network/algorithm/full.hpp>
#include <glas3/array/array_network/algorithm/norm.hpp>

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
    glas3::dense_vector<std::ptrdiff_t> s2( {5, 4, 3} ) ;
    glas3::dense_vector<std::ptrdiff_t> r2( {4, 3} ) ;
    glas3::tensor_train<glas3::dense_array<double>> t2(0, s2, r2) ;
    glas3::randomize( t2, g ) ;
    auto t2_full = glas3::full( t2 ) ;

    std::cout << "Tucker tensor construction (deep copy) from value, shape and rank -> boost::make_shared_noinit apparently returns default initialized??" << '\n' ;
    std::cout << "glas3::Tucker_tensor<glas3::dense_array<double>> t3(1, {4, 3, 5}, {4, 3, 2}) ;" << '\n' << '\n' ;
    glas3::dense_vector<std::ptrdiff_t> s3( {4, 3, 5} ) ;
    glas3::dense_vector<std::ptrdiff_t> r3( {4, 3, 2} ) ;
    glas3::Tucker_tensor<glas3::dense_array<double>> t3(0, s3, r3) ;
    glas3::randomize( t3, g ) ;
    auto t3_full = glas3::full( t3 ) ;

	std::cout << "ttt operation (shallow copy): inner product" << '\n' ;
	std::cout << "auto t1_n = glas3::norm( t1, 2 ) ;" << '\n' ;
	std::cout << "auto t2_n = glas3::norm_2( t2 ) ;" << '\n' ;
	std::cout << "auto t3_n = glas3::norm_2( t3 ) ;" << '\n' ;
	auto t1_n = glas3::norm( t1, 2 ) ;
	auto t2_n = glas3::norm_2( t2 ) ;
	auto t3_n = glas3::norm_2( t3 ) ;
	auto t1_n_ref = glas3::norm_2( t1_full ) ;
	auto t2_n_ref = glas3::norm_2( t2_full ) ;
	auto t3_n_ref = glas3::norm_2( t3_full ) ;

	std::cout << "std::abs( t1_n - t1_n_ref ): " << std::abs( t1_n - t1_n_ref ) << '\n' ;
	std::cout << "std::abs( t2_n - t2_n_ref ): " << std::abs( t2_n - t2_n_ref ) << '\n' ;
	std::cout << "std::abs( t3_n - t3_n_ref ): " << std::abs( t3_n - t3_n_ref ) << '\n' ;
    if ( std::abs( t1_n - t1_n_ref ) > 1e-14 ) return 1 ;
    if ( std::abs( t2_n - t2_n_ref ) > 1e-14 ) return 1 ;
    if ( std::abs( t3_n - t3_n_ref ) > 1e-14 ) return 1 ;

	std::cout << '\n' ;

	return 0;

}
