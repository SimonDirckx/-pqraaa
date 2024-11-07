#include <glas3/array/array_network/type/Tucker_tensor.hpp>

#include <glas3/array/array_network/algorithm/randomize.hpp>
#include <glas3/array/array_network/algorithm/iostream.hpp>

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

 	std::ptrdiff_t i ;

    std::cout << "Tucker tensor construction (deep copy) from value, shape and rank -> boost::make_shared_noinit apparently returns default initialized??" << '\n' ;
    std::cout << "glas3::Tucker_tensor<glas3::dense_array<double>> t(1, {3, 4, 5}, {4, 3, 2}) ;" << '\n' << '\n' ;
    glas3::dense_vector<std::ptrdiff_t> s( {3, 4, 5} ) ;
    glas3::dense_vector<std::ptrdiff_t> r( {4, 3, 2} ) ;
    glas3::Tucker_tensor<glas3::dense_array<double>> t(0, s, r) ;
    glas3::randomize( t, g ) ;

    if ( t.size() != glas3::prod( s ) || t.size() != glas3::prod( s )
    || t.ndof() != glas3::sum( s * r ) + glas3::prod( r )
    || t.ndof() != glas3::sum( s * r ) + glas3::prod( r )
    || t.shape().size() != s.size() || t.shape().shape().size() != 1
    || t.rank().size() != 3
    || ! ( glas3::all( t.rank() == r )
        || glas3::all( t.rank() == r[{0, 2, 1}] )
        || glas3::all( t.rank() == r[{1, 0, 2}] )
        || glas3::all( t.rank() == r[{1, 2, 0}] )
        || glas3::all( t.rank() == r[{2, 0, 1}] )
        || glas3::all( t.rank() == r[{2, 1, 0}] ) ) ) return 1 ;

    for ( i = 0; i < t.shape().size(); ++i ) {
    	if ( t.shape()[i] != s[i] || t.shape()[i] != s[i] ) return 1 ;
    }

    auto t_ref = glas3::full( glas3::ttt_1D( glas3::full( glas3::ttt_1D( glas3::full( glas3::ttt_1D( *t.nodes2arrays().at( 0 ), *t.nodes2arrays().at( 3 ), 2, 1 ) ), *t.nodes2arrays().at( 2 ), 1, 1 ) ), *t.nodes2arrays().at( 1 ), 0, 1 ) ) ;
	glas3::shape_index<std::ptrdiff_t> index( t.shape() ) ;

	for ( index.reset(); index.overflow_count() < 1; index.inc_lin_in() ) {
		if ( std::abs( t_ref( index ) - t( index ) ) > 1e-15
			 || std::abs( t_ref[ index.lin_out() ] - t[ index.lin_out() ] ) > 1e-15
			 || std::abs( t_ref( index ) - t_ref[ index.lin_out() ] ) > 1e-15
			                         ) return 1 ;
	}

	std::cout << t ;

	std::cout << '\n' ;

	return 0;
}
