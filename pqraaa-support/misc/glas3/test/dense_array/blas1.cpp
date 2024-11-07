#include <glas3/bindings/array/dense_array/adaptor_dense_vector.hpp>

#include <glas3/array/dense_array/container/dense_vector.hpp>

#include <glas3/array/dense_array/algorithm/norm.hpp>
#include <glas3/array/dense_array/algorithm/ops.hpp>
#include <glas3/array/dense_array/algorithm/inner_prod.hpp>

#include <boost/numeric/bindings/blas/level1.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	using namespace boost::numeric::bindings ;

 	glas3::dense_vector<double> v({1, 2, 3}), w({4, 5, 6}) ;
 	glas3::dense_vector<double> x = w ;

 	blas::axpy( 2.0, v, w ) ;

 	if (glas3::norm_1( w - glas3::dense_vector<double>(2.0, 3) * v - x ) > 1.e-12 ) return 1 ;

 	if (std::abs( blas::dot( v, w ) - glas3::inner_prod( v, w ) ) > 1.e-12 ) return 1 ;

 	std::cout << '\n' ;

 	return 0 ;

}
