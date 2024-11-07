#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/ops.hpp>

#include <iostream>
#include <cmath>
#include <functional>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j ;

    std::cout << "entrywise_operation" << '\n' ;
    std::cout << "glas3::dense_array<int> a1({1, 2, 3, 4, 5, 6}, {2, 3}) ;" << '\n' ;
    std::cout << "glas3::dense_vector<int> v1({1, 2, 3, 4, 5, 6}) ;" << '\n' ;
    std::cout << "auto b1 = ( a1 * ( v1 + a1 ) - v1 ) / v1 ;" << '\n' ;
    glas3::dense_array<int> a1({1, 2, 3, 4, 5, 6}, {2, 3}) ;
    glas3::dense_vector<int> v1({1, 2, 3, 4, 5, 6}) ;
    auto b1 = ( a1 * ( v1 + a1 ) - v1 ) / v1 ;
    if ( b1.size() != 6 || b1.ndof() != 5 * 6 || b1.shape().size() != 2 || b1.shape()[0] != 2 || b1.shape()[1] != 3 || b1.shape().shape().size() != 1 ) return 1 ;
    a1[0] -= 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    			if ( b1[i + 2 * j] != ( ( a1[i + 2 * j] * ( v1[i + 2 * j] + a1[i + 2 * j] ) - v1[i + 2 * j] ) / v1[i + 2 * j] )
    					|| b1({i, j}) != ( ( a1({i, j}) * ( v1({i + 2 * j}) + a1({i, j}) ) - v1({i + 2 * j}) ) / v1({i + 2 * j}) ) )
    				return 1 ;
    	}
    }

    std::cout << "entrywise_operation" << '\n' ;
    std::cout << "glas3::dense_array<int> a2({1, 2, 3, 4, 5, 6}, {2, 3}) ;" << '\n' ;
    std::cout << "glas3::dense_vector<int> v2({1, 2, 3, 4, 5, 6}) ;" << '\n' ;
    std::cout << "auto b2 = glas3::for_each( std::function< double (double) >( [] (double x) { return std::abs( x ) ; } ), a2 ) ;" << '\n' ;
    std::cout << "auto b3 = glas3::for_each( std::function< double ( double const&, double const& ) >( [] ( double const& x, double const& y ) { return std::pow( x, y ) ; } ), v2, a2 ) ;" << '\n' ;
    glas3::dense_array<int> a2({1, 2, 3, 4, 5, 6}, {2, 3}) ;
    glas3::dense_vector<int> v2({1, 2, 3, 4, 5, 6}) ;
    auto b2 = glas3::for_each( std::function< double (double) >( [] (double x) { return std::abs( x ) ; } ), a2 ) ;
    auto b3 = glas3::for_each( std::function< double ( double const&, double const& ) >( [] ( double const& x, double const& y ) { return std::pow( x, y ) ; } ), v2, a2 ) ;
    if ( b2.size() != 6 || b2.ndof() != 6 || b2.shape().size() != 2 || b2.shape()[0] != 2 || b2.shape()[1] != 3 || b2.shape().shape().size() != 1 ) return 1 ;
    if ( b3.size() != 6 || b3.ndof() != 6 + 6 || b3.shape().size() != 1 || b3.shape()[0] != 6 || b3.shape().shape().size() != 0 ) return 1 ;
    a2[0] += 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		if ( b2[i + 2 * j] != std::abs( a2[i + 2 * j] ) || b2({i, j}) != std::abs( a2({i, j}) ) ) return 1 ;
    		if ( b3[i + 2 * j] != std::pow( v2[i + 2 * j], a2[i + 2 * j] ) || b3({i + 2 * j}) != std::pow( v2({i + 2 * j}), a2({i, j}) ) ) return 1 ;
    	}
    }

    std::cout << "entrywise_operation" << '\n' ;
    std::cout << "glas3::dense_array<double> a3({-3, 0, -2, 0, -1, 9}, {2, 3}) ;" << '\n' ;
    std::cout << "glas3::dense_array<std::complex<double>> a4({-3, 0, std::complex<double>(3, -2), 0, std::complex<double>(0, -1), 9}, {2, 3}) ;" << '\n' ;
    std::cout << "auto a3_conj = glas3::conj( a3 ) ;" << '\n' ;
    std::cout << "auto a4_conj = glas3::conj( a4 ) ;" << '\n' ;
    glas3::dense_array<double> a3({-3, 0, -2, 0, -1, 9}, {2, 3}) ;
    glas3::dense_array<std::complex<double>> a4({-3, 0, std::complex<double>(3, -2), 0, std::complex<double>(0, -1), 9}, {2, 3}) ;
    auto a3_conj = glas3::conj( a3 ) ;
    auto a4_conj = glas3::conj( a4 ) ;
    if ( a3_conj.size() != 6 || a3_conj.ndof() != 6 || a3_conj.shape().size() != 2 || a3_conj.shape()[0] != 2 || a3_conj.shape()[1] != 3 || a3_conj.shape().shape().size() != 1 ) return 1 ;
    if ( a4_conj.size() != 6 || a4_conj.ndof() != 6 || a4_conj.shape().size() != 2 || a4_conj.shape()[0] != 2 || a4_conj.shape()[1] != 3 || a4_conj.shape().shape().size() != 1 ) return 1 ;
    a3[0] += 1 ;
    a4[0] += 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		if ( a3_conj[i + 2 * j] != a3[i + 2 * j] || a3_conj({i, j}) != a3({i, j}) ) return 1 ;
    		if ( a4_conj[i + 2 * j] != std::conj( a4[i + 2 * j] ) || a4_conj({i, j}) != std::conj( a4({i, j}) ) ) return 1 ;
    	}
    }

    std::cout << '\n' ;

	return 0;

}
