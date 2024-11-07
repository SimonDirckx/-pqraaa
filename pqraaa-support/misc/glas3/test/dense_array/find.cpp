#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/find.hpp>
#include <glas3/array/dense_array/algorithm/transpose.hpp>

#include <iostream>

int main() {

	std::cout << '\n';

	std::cout << "find nonzero entries in dense_vector" << '\n' ;
    glas3::dense_vector<double> v = {5, 0, 1, -4, 0, 0} ;
    auto ind = glas3::find( v ) ;
    int c = 0 ;
    for ( int i = 0; i < v.size(); ++i ) {
        if ( v[i] != 0 ) {
        	if ( i != ind[c] ) return 1 ;
        	++c ;
        }
    }
    if ( c != ind.size() ) return 1 ;

    for ( int i = 0; i < ind.size(); ++i ) {
    	std::cout << ind[i] << ' ' ;
    }

    std::cout << "find nonzero entries in (temporary) transposed dense_array" << '\n' ;
    glas3::dense_array<double> a( {5, 0, 1, -4, 0, 0}, {2, 3} ) ;
    auto a_t = glas3::transpose( a ) ;
    ind = glas3::find( glas3::transpose( a ) ) ;
    c = 0 ;
    for ( int i = 0; i < a_t.size(); ++i ) {
        if ( a_t[i] != 0 ) {
        	if ( i != ind[c] ) return 1 ;
        	++c ;
        }
    }
    if ( c != ind.size() ) return 1 ;

    for ( int i = 0; i < ind.size(); ++i ) {
    	std::cout << ind[i] << ' ' ;
    }

    std::cout << '\n' ;

    return 0;
}
