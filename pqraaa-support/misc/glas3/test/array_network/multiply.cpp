#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>
#include <glas3/array/dense_array/type/eye.hpp>

#include <glas3/array/array_network/algorithm/randomize.hpp>
#include <glas3/array/array_network/algorithm/full.hpp>
#include <glas3/array/array_network/algorithm/multiply.hpp>
#include <glas3/array/array_network/algorithm/iostream.hpp>

#include <glas3/dense_array.hpp>

#include <glas2/sparse.hpp>

#include <random>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <iostream>

int main() {

 	std::cout << '\n';

 	std::cout << "multiply operation" << '\n' ;

 	typedef typename boost::uniform_real<double>                        distribution_type ;
 	typedef typename boost::lagged_fibonacci607                         engine_type ;
 	typedef boost::variate_generator< engine_type, distribution_type >  generator_type ;
 	std::random_device rd ;
 	generator_type g( engine_type( rd() ), distribution_type( -1, 1 ) ) ;

 	glas3::dense_vector<double> v1( glas3::no_init(), 4 ), v2( glas3::no_init(), 4 ), v3( glas3::no_init(), 5 ) ;
 	glas3::dense_matrix<double> m1( glas3::no_init(), {5, 4} ), m2( glas3::no_init(), {4, 5} ), m3( glas3::no_init(), {4, 4} ) ;
 	glas3::dense_array<double> a1( glas3::no_init(), {5, 3, 5, 4} ), a2( glas3::no_init(), {4, 3, 4, 5} ) ;
 	glas3::randomize( v1, g ) ;
 	glas3::randomize( v2, g ) ;
 	glas3::randomize( v3, g ) ;
 	glas3::randomize( m1, g ) ;
 	glas3::randomize( m2, g ) ;
 	glas3::randomize( m3, g ) ;
 	glas3::randomize( a1, g ) ;
 	glas3::randomize( a2, g ) ;
 	glas2::coo<double, short int> s( 5, 4, 3 ) ;
 	s.push_back( 0, 0, 1.0 ) ;
 	s.push_back( 1, 2, 2.0 ) ;
 	s.push_back( 3, 0, 3.0 ) ;
 	s.push_back( 4, 2, 4.0 ) ;
 	push_back( s, glas2::speye( 5, 4 ) ) ;
 	push_back( s, -3.5 * glas2::speye( 4, 4 ), glas2::range( 1, 5 ), glas2::all() ) ;
 	glas3::dense_matrix<double> s_full( glas3::no_init(), {5, 4} ) ;
 	s_full({0, 0}) = 1 ;
 	s_full({1, 2}) = 2 ;
 	s_full({3, 0}) = 3 ;
 	s_full({4, 2}) = 4 ;
 	s_full = s_full + glas3::eye<double>( {5, 4} ) ;
 	s_full({1, 0}) -= 3.5 ;
 	s_full({2, 1}) -= 3.5 ;
 	s_full({3, 2}) -= 3.5 ;
 	s_full({4, 3}) -= 3.5 ;

    std::cout << "CP tensor construction (deep copy) from value, shape and rank" << '\n' ;
    std::cout << "glas3::CP_tensor<glas3::dense_array<double>> t1( 0, {3, 4, 5}, 3 ) ;" << '\n' << '\n' ;
    glas3::dense_vector<std::ptrdiff_t> s1( {3, 4, 5} ) ;
    std::ptrdiff_t r1( 3 ) ;
    glas3::CP_tensor<glas3::dense_array<double>> t1( 0, s1, r1 ) ;
    glas3::randomize( t1, g ) ;
    auto t1_full = glas3::full( t1 ) ;

    std::cout << "tensor_train construction (deep copy) from value, shape and rank" << '\n' ;
    std::cout << "glas3::tensor_train<glas3::dense_array<double>> t2(0, {5, 4, 3}, {4, 3}) ;" << '\n' << '\n' ;
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

 	auto t2_m1_11 = glas3::multiply( t2, m1, 1, 1 ) ;
 	auto t2_m1_11_full = glas3::full( t2_m1_11 ) ;
 	auto t2_m1_11_ref = glas3::multiply( t2_full, m1, 1, 1 ) ;

 	auto t2_s_11 = glas3::multiply( t2, s, 1, 1 ) ;
 	auto t2_s_11_full = glas3::full( t2_s_11 ) ;
 	auto t2_s_11_ref = glas3::multiply( t2_full, s, 1, 1 ) ;

 	auto t3_m3_00 = glas3::multiply( t3, m3, 0, 0 ) ;
 	auto t3_m3_00_full = glas3::full( t3_m3_00 ) ;
 	auto t3_m3_00_ref = glas3::multiply( t3_full, m3, 0, 0 ) ;

 	auto t1_m3_10 = glas3::multiply( t1, m3, 1, 0 ) ;
 	auto t1_m3_10_full = glas3::full( t1_m3_10 ) ;
 	auto t1_m3_10_ref = glas3::multiply( t1_full, m3, 1, 0 ) ;

 	auto t2_a2_12 = glas3::multiply( t2, a2, 1, 2 ) ;
 	auto t2_a2_12_full = glas3::full( t2_a2_12 ) ;
 	auto t2_a2_12_ref = glas3::multiply( t2_full, a2, 1, 2 ) ;

 	if ( ! glas3::all( t2_m1_11.shape() == t2_m1_11_ref.shape() ) || glas3::norm_inf( t2_m1_11_full - t2_m1_11_ref ) > 1e-14 ) return 1 ;
 	if ( ! glas3::all( t2_s_11.shape() == t2_s_11_ref.shape() ) || glas3::norm_inf( t2_s_11_full - t2_s_11_ref ) > 1e-14 ) return 1 ;
 	if ( ! glas3::all( t3_m3_00.shape() == t3_m3_00_ref.shape() ) || glas3::norm_inf( t3_m3_00_full - t3_m3_00_ref ) > 1e-14 ) return 1 ;
 	if ( ! glas3::all( t1_m3_10.shape() == t1_m3_10_ref.shape() ) || glas3::norm_inf( t1_m3_10_full - t1_m3_10_ref ) > 1e-14 ) return 1 ;
 	if ( ! glas3::all( t2_a2_12.shape() == t2_a2_12_ref.shape() ) || glas3::norm_inf( t2_a2_12_full - t2_a2_12_ref ) > 1e-14 ) return 1 ;

	std::cout << '\n' ;

	return 0;

}
