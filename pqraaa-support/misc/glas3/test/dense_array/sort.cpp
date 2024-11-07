#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/sort.hpp>
#include <glas3/array/dense_array/algorithm/transpose.hpp>

#include <iostream>

int main() {

	std::cout << '\n';

	std::cout << "sort dense_vector" << '\n' ;
    glas3::dense_vector<double> v = {5, 3, 1, 4, 7, 6} ;
    glas3::sort( v ) ;
    for ( int i = 1; i < v.size(); ++i ) {
        if ( v[i - 1] > v[i] ) return 1 ;
    }

    std::cout << "sort (temporary) transposed dense_array" << '\n' ;
    glas3::dense_array<double> a( {5, 3, 1, 4, 7, 6}, {2, 3} ) ;
    auto a_t = glas3::transpose( a ) ;
    glas3::sort( glas3::transpose( a ) ) ;
    for ( int i = 1; i < a_t.size(); ++i ) {
        if ( a_t[i - 1] > a_t[i] ) return 1 ;
    }

    auto comp = std::function< bool ( double const&, double const& ) >( [] ( double const& x, double const& y ) { return x > y ; } ) ;
    glas3::sort( glas3::transpose( a ), comp ) ;
    for ( int i = 1; i < a_t.size(); ++i ) {
    	if ( a_t[i - 1] < a_t[i] ) return 1 ;
    }

    glas3::sort_ascending( glas3::transpose( a ) ) ;
    for ( int i = 1; i < a_t.size(); ++i ) {
    	if ( a_t[i - 1] > a_t[i] ) return 1 ;
    }

    glas3::sort_descending( glas3::transpose( a ) ) ;
    for ( int i = 1; i < a_t.size(); ++i ) {
    	if ( a_t[i - 1] < a_t[i] ) return 1 ;
    }

//    glas3::dense_array<int> a1( {5, 3, 1, 4, 7, 6}, {2, 3} ) ;
//    auto a1_t = a1.shallow_copy(); //glas3::transpose( a1 ) ;
//    glas3::dense_array<int> a2 = a1_t ;
//    glas3::dense_array<int> ind = glas3::range<>( 0, a1_t.size() ) ;
//
//    //std::cout << "Before sort " << v << std::endl ;
//    glas3::sort( a1, ind ) ;
//   // std::cout << "After sort " << v << std::endl ;
//   //std::cout << ind << std::endl ;
//
//    for ( int i = 1; i < a1_t.size(); ++i ) {
//    	if ( a1_t[i - 1] > a1_t[i] ) return 2 ;
//    }
//    for ( int i = 0; i < a1_t.size(); ++i ) {
//    	if ( a2[ind[i]] != a1_t[i] ) return 3 ;
//    }

    std::cout << '\n' ;

    return 0;
}
