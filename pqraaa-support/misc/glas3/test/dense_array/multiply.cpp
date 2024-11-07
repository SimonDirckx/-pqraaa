#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>
#include <glas3/array/dense_array/type/eye.hpp>

#include <glas3/array/dense_array/algorithm/multiply.hpp>
#include <glas3/array/dense_array/algorithm/randomize.hpp>
#include <glas3/array/dense_array/algorithm/all.hpp>
#include <glas3/array/dense_array/algorithm/ttt_1D.hpp>
#include <glas3/array/dense_array/algorithm/ops.hpp>
#include <glas3/array/dense_array/algorithm/all.hpp>

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
 	glas3::dense_matrix<double> m1( glas3::no_init(), {5, 4} ), m2( glas3::no_init(), {4, 5} ) ;
 	glas3::dense_array<double> a1( glas3::no_init(), {5, 3, 5, 4} ), a2( glas3::no_init(), {4, 3, 4, 5} ) ;
 	glas3::randomize( v1, g ) ;
 	glas3::randomize( v2, g ) ;
 	glas3::randomize( v3, g ) ;
 	glas3::randomize( m1, g ) ;
 	glas3::randomize( m2, g ) ;
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

 	auto v1_v2 = glas3::multiply( v1, v2 ) ;
 	auto v3_m1 = glas3::multiply( v3, m1 ) ;
 	auto m1_v2 = glas3::multiply( m1, v2 ) ;
 	auto m1_m2 = glas3::multiply( m1, m2 ) ;
 	auto a1_v1 = glas3::multiply( a1, v1 ) ;
 	auto m2_a1 = glas3::multiply( m2, a1 ) ;
 	auto s_v2 = glas3::multiply( s, v2 ) ;
 	auto v3_s = glas3::multiply( v3, s ) ;
 	auto s_m2 = glas3::multiply( s, m2 ) ;
 	auto m2_s = glas3::multiply( m2, s ) ;
 	auto s_a2 = glas3::multiply( s, a2 ) ;
 	auto a2_s = glas3::multiply( a2, s ) ;
 	auto v1_v2_00 = glas3::multiply( v1, v2, 0, 0 ) ;
 	auto v3_m1_00 = glas3::multiply( v3, m1, 0, 0 ) ;
 	auto m1_v2_10 = glas3::multiply( m1, v2, 1, 0 ) ;
 	auto m1_m2_10 = glas3::multiply( m1, m2, 1, 0 ) ;
 	auto a1_v1_30 = glas3::multiply( a1, v1, 3, 0 ) ;
 	auto m2_a1_10 = glas3::multiply( m2, a1, 1, 0 ) ;
	auto s_v2_10 = glas3::multiply( s, v2, 1, 0 ) ;
	auto s_v3_00 = glas3::multiply( s, v3, 0, 0 ) ;
	auto v2_s_01 = glas3::multiply( v2, s, 0, 1 ) ;
	auto v3_s_00 = glas3::multiply( v3, s, 0, 0 ) ;
	auto s_m2_10 = glas3::multiply( s, m2, 1, 0 ) ;
	auto s_m2_01 = glas3::multiply( s, m2, 0, 1 ) ;
	auto s_m1_00 = glas3::multiply( s, m1, 0, 0 ) ;
	auto s_m1_11 = glas3::multiply( s, m1, 1, 1 ) ;
	auto m2_s_10 = glas3::multiply( m2, s, 1, 0 ) ;
	auto m2_s_01 = glas3::multiply( m2, s, 0, 1 ) ;
	auto m1_s_00 = glas3::multiply( m1, s, 0, 0 ) ;
	auto m1_s_11 = glas3::multiply( m1, s, 1, 1 ) ;
	auto s_a2_10 = glas3::multiply( s, a2, 1, 0 ) ;
	auto s_a2_03 = glas3::multiply( s, a2, 0, 3 ) ;
	auto s_a2_12 = glas3::multiply( s, a2, 1, 2 ) ;
	auto s_a1_02 = glas3::multiply( s, a1, 0, 2 ) ;
	auto a2_s_01 = glas3::multiply( a2, s, 0, 1 ) ;
	auto a2_s_30 = glas3::multiply( a2, s, 3, 0 ) ;
	auto a2_s_21 = glas3::multiply( a2, s, 2, 1 ) ;
	auto a1_s_20 = glas3::multiply( a1, s, 2, 0 ) ;

 	if ( v1_v2.shape().size() != glas3::ttt_1D( v1, v2 ).shape().size() || ! glas3::all( v1_v2 == glas3::ttt_1D( v1, v2 ) ) ) return 1 ;
    if ( v3_m1.shape().size() != glas3::ttt_1D( v3, m1 ).shape().size() || ! glas3::all( v3_m1.shape() == glas3::ttt_1D( v3, m1 ).shape() ) || ! glas3::all( v3_m1 == glas3::ttt_1D( v3, m1 ) ) ) return 1 ;
 	if ( m1_v2.shape().size() != glas3::ttt_1D( m1, v2 ).shape().size() || ! glas3::all( m1_v2.shape() == glas3::ttt_1D( m1, v2 ).shape() ) || ! glas3::all( m1_v2 == glas3::ttt_1D( m1, v2 ) ) ) return 1 ;
 	if ( m1_m2.shape().size() != glas3::ttt_1D( m1, m2 ).shape().size() || ! glas3::all( m1_m2.shape() == glas3::ttt_1D( m1, m2 ).shape() ) || ! glas3::all( m1_m2 == glas3::ttt_1D( m1, m2 ) ) ) return 1 ;
 	if ( a1_v1.shape().size() != glas3::ttt_1D( a1, v1 ).shape().size() || ! glas3::all( a1_v1.shape() == glas3::ttt_1D( a1, v1 ).shape() ) || ! glas3::all( a1_v1 == glas3::ttt_1D( a1, v1 ) ) ) return 1 ;
    if ( m2_a1.shape().size() != glas3::ttt_1D( m2, a1 ).shape().size() || ! glas3::all( m2_a1.shape() == glas3::ttt_1D( m2, a1 ).shape() ) || ! glas3::all( m2_a1 == glas3::ttt_1D( m2, a1 ) ) ) return 1 ;
 	if ( s_v2.shape().size() != glas3::ttt_1D( s_full, v2 ).shape().size() || ! glas3::all( s_v2.shape() == glas3::ttt_1D( s_full, v2 ).shape() ) || ! glas3::all( s_v2 == glas3::ttt_1D( s_full, v2 ) ) ) return 1 ;
 	if ( v3_s.shape().size() != glas3::ttt_1D( v3, s_full ).shape().size() || ! glas3::all( v3_s.shape() == glas3::ttt_1D( v3, s_full ).shape() ) || ! glas3::all( v3_s == glas3::ttt_1D( v3, s_full ) ) ) return 1 ;
 	if ( s_m2.shape().size() != glas3::ttt_1D( s_full, m2 ).shape().size() || ! glas3::all( s_m2.shape() == glas3::ttt_1D( s_full, m2 ).shape() ) || ! glas3::all( s_m2 == glas3::ttt_1D( s_full, m2 ) ) ) return 1 ;
 	if ( m2_s.shape().size() != glas3::ttt_1D( m2, s_full ).shape().size() || ! glas3::all( m2_s.shape() == glas3::ttt_1D( m2, s_full ).shape() ) || ! glas3::all( m2_s == glas3::ttt_1D( m2, s_full ) ) ) return 1 ;
 	if ( s_a2.shape().size() != glas3::ttt_1D( s_full, a2 ).shape().size() || ! glas3::all( s_a2.shape() == glas3::ttt_1D( s_full, a2 ).shape() ) || ! glas3::all( s_a2 == glas3::ttt_1D( s_full, a2 ) ) ) return 1 ;
 	if ( a2_s.shape().size() != glas3::ttt_1D( a2, s_full ).shape().size() || ! glas3::all( a2_s.shape() == glas3::ttt_1D( a2, s_full ).shape() ) || ! glas3::all( a2_s == glas3::ttt_1D( a2, s_full ) ) ) return 1 ;
 	if ( v1_v2_00.shape().size() != glas3::ttt_1D( v1, v2, 0, 0 ).shape().size() || ! glas3::all( v1_v2_00 == glas3::ttt_1D( v1, v2, 0, 0 ) ) ) return 1 ;
    if ( v3_m1_00.shape().size() != glas3::ttt_1D( v3, m1, 0, 0 ).shape().size() || ! glas3::all( v3_m1_00.shape() == glas3::ttt_1D( v3, m1, 0, 0 ).shape() ) || ! glas3::all( v3_m1_00 == glas3::ttt_1D( v3, m1, 0, 0 ) ) ) return 1 ;
 	if ( m1_v2_10.shape().size() != glas3::ttt_1D( m1, v2, 1, 0 ).shape().size() || ! glas3::all( m1_v2_10.shape() == glas3::ttt_1D( m1, v2, 1, 0 ).shape() ) || ! glas3::all( m1_v2_10 == glas3::ttt_1D( m1, v2, 1, 0 ) ) ) return 1 ;
 	if ( m1_m2_10.shape().size() != glas3::ttt_1D( m1, m2, 1, 0 ).shape().size() || ! glas3::all( m1_m2_10.shape() == glas3::ttt_1D( m1, m2, 1, 0 ).shape() ) || ! glas3::all( m1_m2_10 == glas3::ttt_1D( m1, m2, 1, 0 ) ) ) return 1 ;
 	if ( a1_v1_30.shape().size() != glas3::ttt_1D( a1, v1, 3, 0 ).shape().size() || ! glas3::all( a1_v1_30.shape() == glas3::ttt_1D( a1, v1, 3, 0 ).shape() ) || ! glas3::all( a1_v1_30 == glas3::ttt_1D( a1, v1, 3, 0 ) ) ) return 1 ;
 	if ( m2_a1_10.shape().size() != glas3::ttt_1D( m2, a1, 1, 0 ).shape().size() || ! glas3::all( m2_a1_10.shape() == glas3::ttt_1D( m2, a1, 1, 0 ).shape() ) || ! glas3::all( m2_a1_10 == glas3::ttt_1D( m2, a1, 1, 0 ) ) ) return 1 ;
	if ( s_v2_10.shape().size() != glas3::ttt_1D( s_full, v2, 1, 0 ).shape().size() || ! glas3::all( s_v2_10.shape() == glas3::ttt_1D( s_full, v2, 1, 0 ).shape() ) || ! glas3::all( s_v2_10 == glas3::ttt_1D( s_full, v2, 1, 0 ) ) ) return 1 ;
	if ( s_v3_00.shape().size() != glas3::ttt_1D( s_full, v3, 0, 0 ).shape().size() || ! glas3::all( s_v3_00.shape() == glas3::ttt_1D( s_full, v3, 0, 0 ).shape() ) || ! glas3::all( s_v3_00 == glas3::ttt_1D( s_full, v3, 0, 0 ) ) ) return 1 ;
	if ( v2_s_01.shape().size() != glas3::ttt_1D( v2, s_full, 0, 1 ).shape().size() || ! glas3::all( v2_s_01.shape() == glas3::ttt_1D( v2, s_full, 0, 1 ).shape() ) || ! glas3::all( v2_s_01 == glas3::ttt_1D( v2, s_full, 0, 1 ) ) ) return 1 ;
	if ( v3_s_00.shape().size() != glas3::ttt_1D( v3, s_full, 0, 0 ).shape().size() || ! glas3::all( v3_s_00.shape() == glas3::ttt_1D( v3, s_full, 0, 0 ).shape() ) || ! glas3::all( v3_s_00 == glas3::ttt_1D( v3, s_full, 0, 0 ) ) ) return 1 ;
	if ( s_m2_10.shape().size() != glas3::ttt_1D( s_full, m2, 1, 0 ).shape().size() || ! glas3::all( s_m2_10.shape() == glas3::ttt_1D( s_full, m2, 1, 0 ).shape() ) || ! glas3::all( s_m2_10 == glas3::ttt_1D( s_full, m2, 1, 0 ) ) ) return 1 ;
	if ( s_m2_01.shape().size() != glas3::ttt_1D( s_full, m2, 0, 1 ).shape().size() || ! glas3::all( s_m2_01.shape() == glas3::ttt_1D( s_full, m2, 0, 1 ).shape() ) || ! glas3::all( s_m2_01 == glas3::ttt_1D( s_full, m2, 0, 1 ) ) ) return 1 ;
	if ( s_m1_00.shape().size() != glas3::ttt_1D( s_full, m1, 0, 0 ).shape().size() || ! glas3::all( s_m1_00.shape() == glas3::ttt_1D( s_full, m1, 0, 0 ).shape() ) || ! glas3::all( s_m1_00 == glas3::ttt_1D( s_full, m1, 0, 0 ) ) ) return 1 ;
	if ( s_m1_11.shape().size() != glas3::ttt_1D( s_full, m1, 1, 1 ).shape().size() || ! glas3::all( s_m1_11.shape() == glas3::ttt_1D( s_full, m1, 1, 1 ).shape() ) || ! glas3::all( s_m1_11 == glas3::ttt_1D( s_full, m1, 1, 1 ) ) ) return 1 ;
	if ( m2_s_10.shape().size() != glas3::ttt_1D( m2, s_full, 1, 0 ).shape().size() || ! glas3::all( m2_s_10.shape() == glas3::ttt_1D( m2, s_full, 1, 0 ).shape() ) || ! glas3::all( m2_s_10 == glas3::ttt_1D( m2, s_full, 1, 0 ) ) ) return 1 ;
	if ( m2_s_01.shape().size() != glas3::ttt_1D( m2, s_full, 0, 1 ).shape().size() || ! glas3::all( m2_s_01.shape() == glas3::ttt_1D( m2, s_full, 0, 1 ).shape() ) || ! glas3::all( m2_s_01 == glas3::ttt_1D( m2, s_full, 0, 1 ) ) ) return 1 ;
	if ( m1_s_00.shape().size() != glas3::ttt_1D( m1, s_full, 0, 0 ).shape().size() || ! glas3::all( m1_s_00.shape() == glas3::ttt_1D( m1, s_full, 0, 0 ).shape() ) || ! glas3::all( m1_s_00 == glas3::ttt_1D( m1, s_full, 0, 0 ) ) ) return 1 ;
	if ( m1_s_11.shape().size() != glas3::ttt_1D( m1, s_full, 1, 1 ).shape().size() || ! glas3::all( m1_s_11.shape() == glas3::ttt_1D( m1, s_full, 1, 1 ).shape() ) || ! glas3::all( m1_s_11 == glas3::ttt_1D( m1, s_full, 1, 1 ) ) ) return 1 ;
	if ( s_a2_10.shape().size() != glas3::ttt_1D( s_full, a2, 1, 0 ).shape().size() || ! glas3::all( s_a2_10.shape() == glas3::ttt_1D( s_full, a2, 1, 0 ).shape() ) || ! glas3::all( s_a2_10 == glas3::ttt_1D( s_full, a2, 1, 0 ) ) ) return 1 ;
	if ( s_a2_03.shape().size() != glas3::ttt_1D( s_full, a2, 0, 3 ).shape().size() || ! glas3::all( s_a2_03.shape() == glas3::ttt_1D( s_full, a2, 0, 3 ).shape() ) || ! glas3::all( s_a2_03 == glas3::ttt_1D( s_full, a2, 0, 3 ) ) ) return 1 ;
	if ( s_a2_12.shape().size() != glas3::ttt_1D( s_full, a2, 1, 2 ).shape().size() || ! glas3::all( s_a2_12.shape() == glas3::ttt_1D( s_full, a2, 1, 2 ).shape() ) || ! glas3::all( s_a2_12 == glas3::ttt_1D( s_full, a2, 1, 2 ) ) ) return 1 ;
	if ( s_a1_02.shape().size() != glas3::ttt_1D( s_full, a1, 0, 2 ).shape().size() || ! glas3::all( s_a1_02.shape() == glas3::ttt_1D( s_full, a1, 0, 2 ).shape() ) || ! glas3::all( s_a1_02 == glas3::ttt_1D( s_full, a1, 0, 2 ) ) ) return 1 ;
	if ( a2_s_01.shape().size() != glas3::ttt_1D( a2, s_full, 0, 1 ).shape().size() || ! glas3::all( a2_s_01.shape() == glas3::ttt_1D( a2, s_full, 0, 1 ).shape() ) || ! glas3::all( a2_s_01 == glas3::ttt_1D( a2, s_full, 0, 1 ) ) ) return 1 ;
	if ( a2_s_30.shape().size() != glas3::ttt_1D( a2, s_full, 3, 0 ).shape().size() || ! glas3::all( a2_s_30.shape() == glas3::ttt_1D( a2, s_full, 3, 0 ).shape() ) || ! glas3::all( a2_s_30 == glas3::ttt_1D( a2, s_full, 3, 0 ) ) ) return 1 ;
	if ( a2_s_21.shape().size() != glas3::ttt_1D( a2, s_full, 2, 1 ).shape().size() || ! glas3::all( a2_s_21.shape() == glas3::ttt_1D( a2, s_full, 2, 1 ).shape() ) || ! glas3::all( a2_s_21 == glas3::ttt_1D( a2, s_full, 2, 1 ) ) ) return 1 ;
	if ( a1_s_20.shape().size() != glas3::ttt_1D( a1, s_full, 2, 0 ).shape().size() || ! glas3::all( a1_s_20.shape() == glas3::ttt_1D( a1, s_full, 2, 0 ).shape() ) || ! glas3::all( a1_s_20 == glas3::ttt_1D( a1, s_full, 2, 0 ) ) ) return 1 ;

	std::cout << '\n' ;

	return 0;

}
