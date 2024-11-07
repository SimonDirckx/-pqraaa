#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/norm.hpp>

#include <iostream>
#include <cmath>
#include <complex>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i ;

 	std::cout << "norm operation" << '\n' ;

 	glas3::dense_array<double> a({-3, 0, -2, 0, -1, 9}, {2, 3, 1}) ;
 	glas3::dense_array<std::complex<double>> a2({-3, 0, std::complex<double>(0, -2), 0, std::complex<double>(0, -1), 9}, {2, 3, 1}) ;

 	double dum = 0, dum2 = 0 ;
 	std::complex<double> complex_0(0, 0) ;
    for ( i = 0; i < a.size(); ++i ) {
    	 if ( a[i] != 0 ) dum += 1 ;
    	 if ( a2[i] != complex_0 ) dum2 += 1 ;
    }
    if ( dum != glas3::norm_0(a) || dum != glas3::norm(a, 0) ) return 1 ;
    if ( dum2 != glas3::norm_0(a2) || dum2 != glas3::norm(a2, 0) ) return 1 ;

    dum = 0 ; dum2 = 0 ;
    for ( i = 0; i < a.size(); ++i ) {
    	dum += std::abs( a[i] ) ;
    	dum2 += std::abs( a2[i] ) ;
    }
    if ( dum != glas3::norm_1(a) || dum != glas3::norm(a, 1) ) return 1 ;
    if ( dum2 != glas3::norm_1(a2) || dum2 != glas3::norm(a2, 1) ) return 1 ;

    dum = 0 ; dum2 = 0 ;
    for ( i = 0; i < a.size(); ++i ) {
    	dum += std::pow( a[i], 2 ) ;
    	dum2 += std::norm( a2[i] ) ;
    }
    dum = std::sqrt( dum ) ;
    dum2 = std::sqrt( dum2 ) ;
    if ( dum != glas3::norm_2(a) || dum != glas3::norm(a, 2) ) return 1 ;
    if ( dum2 != glas3::norm_2(a2) || dum2 != glas3::norm(a2, 2) ) return 1 ;

    dum = 0 ; dum2 = 0 ;
    for ( i = 0; i < a.size(); ++i ) {
    	int v( std::abs( a[i] ) ) ;
    	if (v > dum) dum = v ;
    	int v2( std::abs( a2[i] ) ) ;
    	if (v2 > dum2) dum2 = v2 ;
    }
    if ( dum != glas3::norm_inf(a) || dum != glas3::norm(a, INFINITY) ) return 1 ;
    if ( dum2 != glas3::norm_inf(a2) || dum2 != glas3::norm(a2, INFINITY) ) return 1 ;

    dum = 0 ; dum2 = 0 ;
    for ( i = 0; i < a.size(); ++i ) {
    	dum += std::pow( std::abs( a[i] ), 3 ) ;
    	dum2 += std::pow( std::abs( a2[i] ), 3 ) ;
    }
    dum = std::pow ( dum,  1. / 3 ) ;
    dum2 = std::pow ( dum2,  1. / 3 ) ;
    if ( dum != glas3::norm(a, 3) ) return 1 ;
    if ( dum2 != glas3::norm(a2, 3) ) return 1 ;

    std::cout << '\n' ;

	return 0;

}
