#include <glas3/array/array_network/type/CP_tensor.hpp>

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

 	std::ptrdiff_t i, j, k ;

    std::cout << "CP tensor construction (deep copy) from value, shape and rank" << '\n' ;
    std::cout << "glas3::CP_tensor<glas3::dense_array<double>> t( 0, {3, 4, 5}, 3 ) ;" << '\n' << '\n' ;
    glas3::dense_vector<std::ptrdiff_t> s( {3, 4, 5} ) ;
    std::ptrdiff_t r( 3 ) ;
    glas3::CP_tensor<glas3::dense_array<double>> t( 0, s, r ) ;
    glas3::randomize( t, g ) ;

    if ( t.size() != glas3::prod( s ) || t.size() != glas3::prod( s )
    || t.ndof() != glas3::sum( glas3::constant_array<std::ptrdiff_t>( r, s.shape() ) * s )
    || t.ndof() != glas3::sum( glas3::constant_array<std::ptrdiff_t>( r, s.shape() ) * s )
    || t.shape().size() != s.size() || t.shape().shape().size() != 1
    || t.rank().size() != 1 || t.rank()[0] != r ) return 1 ;

    for ( i = 0; i < t.shape().size(); ++i ) {
    	if ( t.shape()[i] != s[i] || t.shape()[i] != s[i] ) return 1 ;
    }

	glas3::shape_index<std::ptrdiff_t> index( t.shape() ) ;

	for ( index.reset(); index.overflow_count() < 1; index.inc_lin_in() ) {
		double sum = 0 ;
		for ( j = 0; j < t.rank()[0]; ++j ) {
			double product = 1 ;
			k = 0 ;
			for ( auto edge: t.outer_edges() ) {
				auto node = std::get<0>( *t.edges2nodes().at(edge).begin() ) ;
				product *= (*t.nodes2arrays().at(node))( { index[k], j } ) ;
				++k ;
			}
			sum += product ;
		}
		if ( sum != t( index ) || sum != t[ index.lin_out() ] ) return 1 ;
	}

	std::cout << t ;

	std::cout << '\n' ;

	return 0;
}
