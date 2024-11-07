#include <glas3/glas3.hpp>

#include <boost/numeric/bindings/lapack/driver/gesv.hpp>

#include <iostream>
#include <cmath>

int main() {

 	std::cout << '\n';

 	glas3::dense_array<double> a(glas3::no_init(), {2, 3, 2}) ;
 	a = glas3::range<>(0, 12) ;
	glas3::dense_vector<int> v = {1, 2, 0, 1, 5} ;
 	auto a_v = a[v] ;

 	std::cout << "glas3::dense_array<int> a({0, ..., 11}, {2, 3, 2}) ;" << '\n' ;
 	std::cout << "glas3::dense_vector<int> v = {1, 2, 0, 1, 5} ;" << '\n' << '\n' ;

 	std::cout << "auto a_v = a[v] ;" << '\n' << a[v] << '\n' ;

 	std::cout << "glas3::dense_vector<int> a_v_deep_copy = a[v] ;" << '\n' << '\n' ;

 	a_v[2] = -1 ;
 	std::cout << "a_v[2] = -1 ;" << '\n' ;
 	std::cout << "a ;" << '\n' << a << '\n' ;

    auto a_m = a({0, glas3::range<>(1, 3), {1, 0}}) ;
    std::cout << "auto a_m = a({0, glas3::range<>(1, 3), {1, 0}}) ;" << '\n' << a_m << '\n' ;

    std::cout << "other view operations:" << '\n' << "[i]permute, transpose, reshape, vect, squeeze, tile, ttt, ttt_1D, +, -, *, /, &&, <, !=, abs, conj, cos, arg, log, exp, ..." << '\n' << '\n' ;

    std::cout << "some examples:" << '\n' << '\n' ;

    glas3::sort( glas3::transpose( a_m ) ) ;
    std::cout << "glas3::sort( glas3::transpose( a_m ) ) ;" << '\n' << glas3::transpose( a_m ) ;
    std::cout << "a_m ;" << '\n' << a_m ;
    std::cout << "a ;" << '\n' << a << '\n' ;

    glas3::any( a == a * a ) ;
    std::cout << "glas3::any( a == a * a ) ;" << '\n' << glas3::any( a == a * a ) << '\n' << '\n' ;

    glas3::dense_vector<double> X({1, 2, 3}) ;
    glas3::dense_vector<int> ind({7, 8, 9}) ;
    glas3::dense_matrix<double> A({1, 2, 3, 0, 5, 9, -1, 13, 9}, {3, 3}) ;
    auto B = glas3::ttt(A, X, {1}, {0}) ;
    std::cout << "glas3::dense_vector<double> X({1, 2, 3}) ;" << '\n' ;
    std::cout << "glas3::dense_matrix<double> A({1, 2, 3, 0, 5, 9, -1, 13, 9}, {3, 3}) ;" << '\n' ;
    std::cout << "auto B = glas3::ttt(A, X, {1}, {0}) ;" << '\n' << B << '\n' ;

    using namespace boost::numeric::bindings ;
    glas3::dense_vector<double> B2 = B ;
    lapack::gesv( A, ind, B2 ) ;
    std::cout << "glas3::dense_vector<double> B2 = B ;" << '\n' ;
    std::cout << "lapack::gesv( A, ind, B2 ) ;" << '\n' ;
    std::cout << "B2 ;" << '\n' << B2 << '\n' ;

//
// 	std::cout << "glas3::dense_array<std::complex<double>> a2({-3, 0, std::complex<double>(0, -2), 0, std::complex<double>(0, -1), 9}, {3, 2}) ;" << '\n' << a2 << '\n' << '\n' ;
// 	std::cout << "glas3::norm(a1, 0) ;" << '\n' << glas3::norm(a1, 0) << '\n' << '\n' ;
// 	std::cout << "glas3::norm_2(a2) ;" << '\n' << glas3::norm_2(a2) << '\n' << '\n' ;
// 	std::cout << "glas3::norm(a2, INFINITY) ;" << '\n' << glas3::norm(a2, INFINITY) << '\n' << '\n' ;
// 	std::cout << "glas3::inner_prod(a1, a2) ;" << '\n' << glas3::inner_prod(a1, a2) << '\n' << '\n' ;
// 	std::cout << "other operations:" << '\n' << "all, any, assign, fill, ndims, ndof, prod, shape, size, sort, sum" << '\n' << '\n' ;

	return 0;
}
