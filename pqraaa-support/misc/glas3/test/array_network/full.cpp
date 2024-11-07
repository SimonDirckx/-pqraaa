#include <glas3/array/array_network/type/CP_tensor.hpp>
#include <glas3/array/array_network/type/Tucker_tensor.hpp>
#include <glas3/array/array_network/type/tensor_train.hpp>

#include <glas3/array/array_network/algorithm/randomize.hpp>
#include <glas3/array/array_network/algorithm/full.hpp>

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

 	std::ptrdiff_t j, k ;

    std::cout << "CP tensor construction (deep copy) from value, shape and rank" << '\n' ;
    std::cout << "glas3::CP_tensor<glas3::dense_array<double>> t( 0, {3, 4, 5}, 3 ) ;" << '\n' << '\n' ;
    glas3::dense_vector<std::ptrdiff_t> s1( {3, 4, 5} ) ;
    std::ptrdiff_t r1( 3 ) ;
    glas3::CP_tensor<glas3::dense_array<double>> t1( 0, s1, r1 ) ;
    glas3::randomize( t1, g ) ;

	glas3::dense_array<double> t_ref1( glas3::no_init(), t1.shape() ) ;
	glas3::shape_index<std::ptrdiff_t> index( t1.shape() ) ;

	for ( index.reset(); index.overflow_count() < 1; index.inc_lin_in() ) {
		double sum = 0 ;
		for ( j = 0; j < t1.rank()[0]; ++j ) {
			double product = 1 ;
			k = 0 ;
			for ( auto edge: t1.outer_edges() ) {
				auto node = std::get<0>( *t1.edges2nodes().at(edge).begin() ) ;
				product *= (*t1.nodes2arrays().at(node))( { index[k], j } ) ;
				++k ;
			}
			sum += product ;
		}
		t_ref1( index ) = sum ;
	}

	if ( glas3::norm_inf( glas3::full( t1 ) - t_ref1 ) > 1e-15 ) return 1 ;

    std::cout << "tensor_train construction (deep copy) from value, shape and rank" << '\n' ;
    std::cout << "glas3::tensor_train<glas3::dense_array<double>> t(0, {3, 4, 5}, {4, 3}) ;" << '\n' ;
    glas3::dense_vector<std::ptrdiff_t> s2( {3, 4, 5} ) ;
    glas3::dense_vector<std::ptrdiff_t> r2( {4, 3} ) ;
    glas3::tensor_train<glas3::dense_array<double>> t2(0, s2, r2) ;
    glas3::randomize( t2, g ) ;

    auto t_ref2 = glas3::full( glas3::ttt_1D( glas3::full( glas3::ttt_1D( *t2.nodes2arrays().at( 0 ), *t2.nodes2arrays().at( 1 ), 1, 0 ) ), *t2.nodes2arrays().at( 2 ), 2, 0 ) ) ;

    if ( glas3::norm_inf( glas3::full( t2 ) - t_ref2 ) > 1e-15 ) return 1 ;

    std::cout << "Tucker tensor construction (deep copy) from value, shape and rank -> boost::make_shared_noinit apparently returns default initialized??" << '\n' ;
    std::cout << "glas3::Tucker_tensor<glas3::dense_array<double>> t(1, {3, 4, 5}, {4, 3, 2}) ;" << '\n' << '\n' ;
    glas3::dense_vector<std::ptrdiff_t> s3( {3, 4, 5} ) ;
    glas3::dense_vector<std::ptrdiff_t> r3( {4, 3, 2} ) ;
    glas3::Tucker_tensor<glas3::dense_array<double>> t3(0, s3, r3) ;
    glas3::randomize( t3, g ) ;

    auto t_ref3 = glas3::full( glas3::ttt_1D( glas3::full( glas3::ttt_1D( glas3::full( glas3::ttt_1D( *t3.nodes2arrays().at( 0 ), *t3.nodes2arrays().at( 3 ), 2, 1 ) ), *t3.nodes2arrays().at( 2 ), 1, 1 ) ), *t3.nodes2arrays().at( 1 ), 0, 1 ) ) ;

    if ( glas3::norm_inf( glas3::full( t3 ) - t_ref3 ) > 1e-15 ) return 1 ;

	std::cout << '\n' ;

	return 0;
}
